#ifndef _ADT3DSHMNEW_HPP_
#define _ADT3DSHMNEW_HPP_

#include <stack>

class Adt3dShmNew {
public:

  // consider replacing with Leafdata: same size, types, etc...
  // for the user, this is more explicit on the ist member...
  
  typedef struct {
    double bbminmax[6];
    int ist,dummy;
  } Bboxdata;
  
private:

  double global_bbmid[3];
  double global_bbrange_max;
  
  typedef struct {
    double bbminmax[6];
    int index[2];
  } Leafdata;

  MPI_Datatype MPI_Leafdata;
  
  int nbb_serial_check;
  double wtime0;

  int nleaves;
  Leafdata * leafdata;
  
#define T Leafdata
# include "Mmap_Munmap.hpp"
#undef T

public:

  Adt3dShmNew(vector<Bboxdata>& mybbVec) {
  
    if (mpi_rank == 0) cout << "Adt3dShmNew()" << endl;

    const int my_nbb = mybbVec.size();
    
    wtime0 = MPI_Wtime();
  
    nbb_serial_check = -1;
    nleaves = 0;
    leafdata = NULL;
    
    // define the Leafdata struct...
    assert(int_size   == 4ll);
    assert(sizeof(Leafdata) == 56);
    assert(sizeof(Bboxdata) == 56);
    
    MPI_Datatype oldtypes[2];
    int blockcounts[2];
    MPI_Aint offsets[2];
    
    // double bbminmax[6]...
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;
    blockcounts[0] = 6;

    // int index[2]...
    offsets[1] = offsets[0] + double_size*6;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 2;
    
    // build the Leafdata type...
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &MPI_Leafdata);
    MPI_Type_commit(&MPI_Leafdata);

    // check that the type size matches the struct size...
    int leafdata_size;
    MPI_Type_size(MPI_Leafdata,&leafdata_size);
    assert(leafdata_size == sizeof(Leafdata));
    
    // Compute the global bounding box...
    double my_bbminmax[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
    for (int ibb = 0; ibb < my_nbb; ++ibb) {
      FOR_I6 my_bbminmax[i] = min(my_bbminmax[i],mybbVec[ibb].bbminmax[i]);
    }
    double bbminmax[6];
    MPI_Allreduce(my_bbminmax,bbminmax,6,MPI_DOUBLE,MPI_MIN,mpi_comm);

    // instead of using a double tol, use an int conversion strategy...

    FOR_I3 global_bbmid[i] = 0.5*(bbminmax[i]-bbminmax[i+3]);
    global_bbrange_max = max(-bbminmax[3]-bbminmax[0],max(-bbminmax[4]-bbminmax[1],-bbminmax[5]-bbminmax[2]));
    assert(global_bbrange_max > 0.0);

    // and the global count...
    int8 my_nbb_int8 = my_nbb;
    int8 nbb_int8;
    MPI_Allreduce(&my_nbb_int8,&nbb_int8,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert(nbb_int8 < TWO_BILLION);
    const int nbb = nbb_int8;
    assert(nbb >= 1);
    
    if (mpi_rank == 0) cout << " > got global nbb: " << nbb << " with bbox: " << 
                         bbminmax[0] << ":" << -bbminmax[3] << " " << 
                         bbminmax[1] << ":" << -bbminmax[4] << " " << 
                         bbminmax[2] << ":" << -bbminmax[5] << endl;

    // ------------------------------------------------------------------
    // figure out the global number of leaves required in the
    // shared memory array on rank 0 of each node (mpi_comm_shared)...
    // ------------------------------------------------------------------
    
    nleaves = 1;
    while (nleaves < nbb) {
      // my_nleaves needs the largest power of 2 greater than n...
      nleaves *= 2;
    }
    // and one final multiplier to allow all leaves to
    // contain just one member...
    nleaves *= 2;
    nleaves -= 1;
    assert(nleaves < TWO_BILLION);
    
    // --------------------------------------------------------------------
    // if this adt is smaller that some threshold size, then 
    // reduce everything to the serial adt build...
    // --------------------------------------------------------------------
    
    // figure out how big our part of the leaf data will be...
    
    int my_nleaves = 0;
    int my_nbb_split = nbb;
    {
      int my_size_split = mpi_size;
      int my_rank_split = mpi_rank;
      while (my_size_split > 1) {
        // for a split size of 1, rank 0 completes the leaf and we are done...
        if (my_nbb_split == 1) {
          if (my_rank_split == 0) {
            my_nleaves += 2;
          }
          my_nbb_split = 0;
          break;
        }
        assert(my_nbb_split > 1);
        // for a mpi_comm size of more than 1, add 2 entries to the leaf data
        // on rank 0 only. The first is for telling the final rank 0 (gathering process)
        // how to unpack, and the second is the actual leafdata...
        if (my_rank_split == 0) {
          my_nleaves += 2;
        }
        int my_size_half = my_size_split/2;
        int my_nbb_half = my_nbb_split/2; 
        if (my_rank_split < my_size_half) {
          my_size_split = my_size_half;
          my_nbb_split = my_nbb_half;
        }
        else {
          my_size_split = my_size_split - my_size_half;
          my_rank_split -= my_size_half;
          my_nbb_split = my_nbb_split - my_nbb_half;
        }
      }
    }
    
    // this check is for the serial function. Get rid of it eventually?
    nbb_serial_check = my_nbb_split;
    
    // on the serial process, use the same counting algo...
    if (my_nbb_split > 0) {
      int my_nleaves_serial = 1;
      while (my_nleaves_serial < my_nbb_split) {
        // my_nleaves needs the largest power of 2 greater than n...
        my_nleaves_serial *= 2;
      }
      // and one final multiplier to allow all leaves to
      // contain just one member...
      my_nleaves_serial *= 2;
      my_nleaves_serial -= 1;
      // add this + 1 to my_nleaves (1 for the indicator)
      my_nleaves += 1 + my_nleaves_serial;
    }
    
    // allocate my_leafdata...
    Leafdata * my_leafdata = new Leafdata[my_nleaves];
    // for now, fill with -2 for checking...
    for (int ileaf = 0; ileaf < my_nleaves; ++ileaf) {
      my_leafdata[ileaf].index[0] = -2;
      my_leafdata[ileaf].index[1] = -2;
    }
    
    // now call the recursive split and shuffle alorithm...
    initRecursive(my_leafdata,mybbVec,my_nbb,nbb,bbminmax,0,mpi_comm);

    // ----------------------------------------------------------------------
    // we (FH,MV) tried another approach to this where we gather each node's
    // leaf data to rank0 of
    // each node, then allgatherx to all rank0's of all nodes, then each
    // does the shuffle in parallel, avoiding the final bcast from global rank0.
    // unfortunates, it was slower when it worked, and failed when large.
    // see code at the bottom of this file.
    // ----------------------------------------------------------------------
    
    // now gather the leaf data at global rank 0...
    int * nleaves_of_rank = NULL;
    if (mpi_rank == 0) nleaves_of_rank = new int[mpi_size];
    MPI_Gather(&my_nleaves,1,MPI_INT,nleaves_of_rank,1,MPI_INT,0,mpi_comm);
    
    int * disp_of_rank = NULL;
    Leafdata * leafdata_buf = NULL;
    int leafdata_count = -1;
    if (mpi_rank == 0) {
      // check for overflow...
      int8 leafdata_count_int8 = 0;
      FOR_RANK leafdata_count_int8 += nleaves_of_rank[rank];
      assert(leafdata_count_int8 < TWO_BILLION);
      disp_of_rank = new int[mpi_size];
      disp_of_rank[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp_of_rank[rank] = disp_of_rank[rank-1] + nleaves_of_rank[rank-1];
      leafdata_count = disp_of_rank[mpi_size-1] + nleaves_of_rank[mpi_size-1];
      leafdata_buf = new Leafdata[leafdata_count];
    }
    
    MPI_Gatherv(my_leafdata,my_nleaves,MPI_Leafdata,leafdata_buf,nleaves_of_rank,disp_of_rank,MPI_Leafdata,0,mpi_comm);
    delete[] my_leafdata;
    
    if (mpi_rank == 0) {
      delete[] nleaves_of_rank;
      delete[] disp_of_rank;
    }
    
    if (mpi_rank == 0) cout << " > done gather: " << MPI_Wtime()-wtime0 << endl;
    
    // everyone allocates the leafdata in shared memory...
    assert(leafdata == NULL);
    CTI_Mmap(leafdata,nleaves);

    if (mpi_rank == 0) cout << " > done shm allocation: " << MPI_Wtime()-wtime0 << endl;
    
    // only rank 0 unpacks into the shared memory. We will bcast later
    // to all nodes...
    if (mpi_rank == 0) {
      
      // for checking...
      //for (int ileaf = 0; ileaf < nleaves; ++ileaf)
      //  leafdata[ileaf].index[0] = -1;

      // unpack the leafdata_buf into the shared memory leafdata...
      int ileafdata = 0;
      while (ileafdata < leafdata_count) {
        
        // the first leaf from any rank should be a "special" leaf 
        // that tells us how to unpack (use -1 indexing)... 
        assert(leafdata_buf[ileafdata].index[0] < 0);
        int ileaf = -leafdata_buf[ileafdata].index[0]-1;
        int nleaf = leafdata_buf[ileafdata].index[1]; // nleaf = 1,3,7,15,...2^x-1
        ++ileafdata;
        
        int nlevel = 1; // 1,2,4,8,16,...
        while (nleaf > 0) {
          // pack nlevel leafdata elements into the ileaf location... 
          for (int ilevel = 0; ilevel < nlevel; ++ilevel) {
            assert(ileaf+ilevel < nleaves);
            //assert(leafdata[ileaf+ilevel].index[0] == -1);
            leafdata[ileaf+ilevel] = leafdata_buf[ileafdata+ilevel];
          }
          // get ready for the next level...
	  ileafdata += nlevel;
          nleaf -= nlevel;
          nlevel *= 2;
          ileaf = 2*ileaf+1;
        }
        // we should have gotten to exactly 0...
        assert(nleaf == 0);
        
      }
      
      delete[] leafdata_buf;
    
    }

    if (mpi_rank == 0) cout << " > done repacking: " << MPI_Wtime()-wtime0 << endl;
    
    // broadcast leafdata from rank0 to all nodes...
    if (mpi_rank_shared == 0) {
      // split bcast...
      const int bcast_size = 2000000000/sizeof(Leafdata);
      int leafdata_offset = 0;
      int count = 0;
      while (leafdata_offset < nleaves) {
        MPI_Bcast(leafdata+leafdata_offset,min(bcast_size,nleaves-leafdata_offset),MPI_Leafdata,0,mpi_comm_internode);
        leafdata_offset += bcast_size;
	++count;
      }
      if (mpi_rank == 0) cout << " > split bcast into " << count << " successive calls" << endl;
      // too big...
      //MPI_Bcast(leafdata,nleaves,MPI_Leafdata,0,mpi_comm_internode);
    }
    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank == 0) cout << " > done bcast: " << MPI_Wtime()-wtime0 << endl;
    
    MPI_Type_free(&MPI_Leafdata);
    
  }

  ~Adt3dShmNew() {
    
    // free the shm...
    CTI_Munmap(leafdata,nleaves);
    
  }
  
  // public check function...
  
  void getRandomIstAndBbox(int& ist,double bbmin[3],double bbmax[3]) {
    
    int ileaf = 0;
    // traverse to the bottom...
    while (leafdata[ileaf].index[0] != leafdata[ileaf].index[1]) {
      // choose a random direction...
      int id = rand()%2; // 0 or 1
      ileaf = 2*ileaf+1+id;
    }
    ist = leafdata[ileaf].index[0];
    FOR_I3 bbmin[i] = leafdata[ileaf].bbminmax[i];
    FOR_I3 bbmax[i] = -leafdata[ileaf].bbminmax[i+3];
    
  }

  void getGlobalBbox(double bbmin[3],double bbmax[3]) const {
    
    // the global bbox is the bbox of the first leaf...
    assert(nleaves > 0);
    FOR_I3 bbmin[i] = leafdata[0].bbminmax[i];
    FOR_I3 bbmax[i] = -leafdata[0].bbminmax[i+3];
    
  }
  
  void addBboxesForPoint(vector<int>& bboxVec,const double x[3]) const {
    
    if (nleaves > 0) {

      stack<int> stack;
      stack.push(0);

      while (!stack.empty()) {
      
        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();

        // skip if out of bounds...
        if ( (x[0] > -leafdata[ileaf].bbminmax[3]) || (x[0] < leafdata[ileaf].bbminmax[0]) ||
             (x[1] > -leafdata[ileaf].bbminmax[4]) || (x[1] < leafdata[ileaf].bbminmax[1]) ||
             (x[2] > -leafdata[ileaf].bbminmax[5]) || (x[2] < leafdata[ileaf].bbminmax[2]) )
          continue;
        
        // this leaf's bbox overlaps. It is the terminal leaf if index is equal
        if (leafdata[ileaf].index[0] == leafdata[ileaf].index[1]) {
          // add this index...
          bboxVec.push_back(leafdata[ileaf].index[0]);
        }
        else {        
          // push both children onto the stack...
          stack.push(2*ileaf+1);
          stack.push(2*ileaf+2);
        }

      }

    }
    
  }

  void addBboxesForSphere(vector<int>& bboxVec,const double x[3],const double r) const {
    
    if (nleaves > 0) {
      
      stack<int> stack;
      stack.push(0);
      
      const double r2 = r*r;
      while (!stack.empty()) {
        
        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();
        
	// calculate the distance of jk to the bbox...
	double d2 = 0;
	FOR_I3 {
	  const double d = max(0.0,max(leafdata[ileaf].bbminmax[i]-x[i],x[i]+leafdata[ileaf].bbminmax[3+i]));
	  d2 += d*d;
	}
        
        // skip if out of bounds...
        if (d2 > r2) 
          continue;
        
        // this leaf's bbox is close enough. It is the terminal leaf if index is equal
        if (leafdata[ileaf].index[0] == leafdata[ileaf].index[1]) {
          // add this index...
          bboxVec.push_back(leafdata[ileaf].index[0]);
        }
        else {        
          // push both children onto the stack...
          stack.push(2*ileaf+1);
          stack.push(2*ileaf+2);
        }

      }

    }
    
  }

  void addBboxesForBbox(vector<int>& bboxVec,const double bbmin[3],const double bbmax[3]) const {
    
    if (nleaves > 0) {
      
      stack<int> stack;
      stack.push(0);
      
      while (!stack.empty()) {

        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();
        
        // do bbox elimination for the leaf...
        if ( (bbmin[0] > -leafdata[ileaf].bbminmax[3]) || (bbmax[0] < leafdata[ileaf].bbminmax[0]) ||
             (bbmin[1] > -leafdata[ileaf].bbminmax[4]) || (bbmax[1] < leafdata[ileaf].bbminmax[1]) ||
             (bbmin[2] > -leafdata[ileaf].bbminmax[5]) || (bbmax[2] < leafdata[ileaf].bbminmax[2]) )
          continue;

        // this leaf's bbox overlaps. It is the terminal leaf if index is equal
        if (leafdata[ileaf].index[0] == leafdata[ileaf].index[1]) {
          // add this index...
          bboxVec.push_back(leafdata[ileaf].index[0]);
        }
        else {        
          // push both children onto the stack...
          stack.push(2*ileaf+1);
          stack.push(2*ileaf+2);
        }

      }

    }
    
  }

  void addBboxesForCylinder(vector<int>& bboxVec,const double x[3],const double n[3],const double r) const {

    if (nleaves > 0) {
      
      stack<int> stack;
      stack.push(0);

      // need to consider a ton of cases here some day...
      // for now, insist on x...
      
      assert(n[0] == 1.0);
      assert(n[1] == 0.0);
      assert(n[2] == 0.0);

      const double r2 = r*r;
      while (!stack.empty()) {
        
        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();
        
	// calculate the distance in y,z of jk to the bbox...
	double d2 = 0;
	for (int i = 1; i < 3; ++i) {
	  const double d = max(0.0,max(leafdata[ileaf].bbminmax[i]-x[i],x[i]+leafdata[ileaf].bbminmax[3+i]));
	  d2 += d*d;
	}
        
        // skip if out of bounds...
        if (d2 > r2) 
          continue;
        
        // this leaf's bbox is close enough. It is the terminal leaf if index is equal
        if (leafdata[ileaf].index[0] == leafdata[ileaf].index[1]) {
          // add this index...
          bboxVec.push_back(leafdata[ileaf].index[0]);
        }
        else {        
          // push both children onto the stack...
          stack.push(2*ileaf+1);
          stack.push(2*ileaf+2);
        }
        
      }
      
    }

  }

  inline double calcMaxD2PointToBbox(const double point[3],const double bbminmax[6]) const {
    double d2 = 0.0;
    FOR_I3 {
      const double dx = max(fabs(bbminmax[i]-point[i]),fabs(-bbminmax[i+3]-point[i]));
      d2 += dx*dx;
    }
    return d2;
  }

  inline double calcMinD2PointToBbox(const double point[3],const double bbminmax[6]) const {
    double d2 = 0.0;
    FOR_I3 {
      if (point[i] < bbminmax[i]) 
	d2 += (bbminmax[i] - point[i])*(bbminmax[i] - point[i]);
      else if (point[i] > -bbminmax[i+3]) 
	d2 += (point[i] + bbminmax[i+3])*(point[i] + bbminmax[i+3]);
    }
    return d2;
  }

  void addBboxesForClosestPoint(vector<int>& bboxVec,const double point[3]) const {
    
    bboxVec.clear();
    
    if (nleaves > 0) {
      
      // each leaf stores its index and d2_min...
      std::list< std::pair<int,double> > leafList;
      leafList.push_back( std::pair<int,double>(0,calcMinD2PointToBbox(point,leafdata[0].bbminmax)) );
      
      // and initialize d2_max with the full bounding box...
      double d2_max = calcMaxD2PointToBbox(point,leafdata[0].bbminmax);
      
      list< std::pair<int,double> >::iterator li = leafList.begin();
      while (li != leafList.end()) {
        
        // recall the d2_min is in second...
        if (li->second <= d2_max) {
          
          // this leaf is valid. If it is a terminal leaf, then increment the iterator...
          if (leafdata[li->first].index[0] == leafdata[li->first].index[1]) {
            
            ++li;
            
          }
          else {

            // otherwise, consider the children...
            const int child1 = 2*(li->first) + 1;
            const int child2 = 2*(li->first) + 2;

            // push the children onto the end of the list -- do not even do the comparison at this point...
            leafList.push_back( std::pair<int,double>(child1,calcMinD2PointToBbox(point,leafdata[child1].bbminmax)) );
            leafList.push_back( std::pair<int,double>(child2,calcMinD2PointToBbox(point,leafdata[child2].bbminmax)) );
            
            // delete the parent...
            list< std::pair<int,double> >::iterator li_copy = li;
            ++li;
            leafList.erase(li_copy);
            
            // and check if the child resets d2_max...
            const double child1_d2_max = calcMaxD2PointToBbox(point,leafdata[child1].bbminmax);
            const double child2_d2_max = calcMaxD2PointToBbox(point,leafdata[child2].bbminmax);
            const double min_d2_max = min(child1_d2_max,child2_d2_max);
            if (min_d2_max < d2_max) {
              
              // the d2_max needs to be changed, and all candidates revisited...
              d2_max = min_d2_max;
              li = leafList.begin();
              
            }
            
          }
          
        }
        else {

          // this leaf now fails, so remove it...
          list< std::pair<int,double> >::iterator li_copy = li;
          ++li;
          leafList.erase(li_copy);
          
        }
       
      }
      
      bboxVec.resize(leafList.size());
      int ii = 0;
      for (li = leafList.begin(); li != leafList.end(); ++li)
        bboxVec[ii++] = leafdata[li->first].index[0];
      
    }
    
  }
  
private:

  void initRecursive(Leafdata *leafdata,vector<Bboxdata>& mybbVec,const int my_nbb,const int nbb,double bbminmax[6],const int ileaf,const MPI_Comm& mpi_split_comm) {
    
    int mpi_split_size;
    MPI_Comm_size(mpi_split_comm,&mpi_split_size);
    
    if (mpi_split_size > 1) {

      int mpi_split_rank;
      MPI_Comm_rank(mpi_split_comm, &mpi_split_rank);

      // ----------------------------------------------------
      // as soon as we are down to 1 bbox amongst ranks, even if
      // there are still a bunch of ranks, stop recursion...
      // ----------------------------------------------------
      
      if (nbb == 1) {
        
        // rank 0 may not have the bbox data, so reduce...
        
        int my_index = -1;
        if (my_nbb == 1) {
          my_index = mybbVec[0].ist;
        }
        else {
          assert(my_nbb == 0);
        }
        int index;
        MPI_Reduce(&my_index,&index,1,MPI_INT,MPI_MAX,0,mpi_split_comm);
        if (mpi_split_rank == 0) {
          
          // make sure we got a valid bbox index...
          assert(index >= 0);
          
          FOR_I2 assert(leafdata[0].index[i] == -2);
          
          // indicator leaf...
          leafdata[0].index[0] = -ileaf-1; // indicates a "special" leaf data...
          leafdata[0].index[1] = 1;
          
          FOR_I2 assert(leafdata[1].index[i] == -2);
          
          leafdata[1].index[0] = leafdata[1].index[1] = index;
          FOR_I6 leafdata[1].bbminmax[i] = bbminmax[i];
          
        }

        // don't go any deeper...
        return;
        
      }

      // ----------------------------------------------------
      // we have more than 1 bbox, so continue to divide 
      // and concur...
      // ----------------------------------------------------
      
      double wtime = MPI_Wtime();
      
      // reduce to get nbb...
      const int nbb_half = nbb/2;
      
      // if we are rank0 on this set of data, we are responsible for adding its
      // info to the leaf data we will be returning to master...
      
      int leafdata_offset = 0;
      if (mpi_split_rank == 0) {
        
        // used below...
        leafdata_offset = 2;
        
        FOR_I2 assert(leafdata[0].index[i] == -2);
        
        // indicator leaf...
        leafdata[0].index[0] = -ileaf-1; // indicates a "special" leaf data...
        leafdata[0].index[1] = 1;
        
        FOR_I2 assert(leafdata[1].index[i] == -2);
        
        // bbox leaf...
        FOR_I6 leafdata[1].bbminmax[i] = bbminmax[i];
        leafdata[1].index[0] = 0; 
        leafdata[1].index[1] = nbb-1; // HACK was my_nbb-1; 
        
        assert(leafdata[1].index[1] != leafdata[1].index[0]);
        
      }
      
      // store the minimum bounds of all boxes as integers...
      
      int my_bbminmax_int[6] = { TWO_BILLION, TWO_BILLION, TWO_BILLION, TWO_BILLION, TWO_BILLION, TWO_BILLION };
      int * buf_int = new int[my_nbb*3];
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        FOR_I3 {
          // store an integer version of the minimum...
          buf_int[3*ibb+i] = (int)floor((mybbVec[ibb].bbminmax[i] - global_bbmid[i])*2.0E+9/global_bbrange_max); // x min, y min, z min as an int
          assert(buf_int[3*ibb+i] < TWO_BILLION);
          assert(buf_int[3*ibb+i] > -TWO_BILLION);
          my_bbminmax_int[i] = min(my_bbminmax_int[i],buf_int[3*ibb+i]);
          my_bbminmax_int[i+3] = min(my_bbminmax_int[i+3],-buf_int[3*ibb+i]);
        }
      }
      int bbminmax_int[6];
      MPI_Allreduce(my_bbminmax_int,bbminmax_int,6,MPI_INT,MPI_MIN,mpi_split_comm);
      
      // because of the rounding associated with the integers, there will be cases when 
      // the upper limit is the correct solution (i.e. contains many of the points), but
      // it cannot be obtained by the bisection routine, so increase this by 1...
      
      FOR_I3 bbminmax_int[i+3] -= 1;
      
      //if (mpi_rank == 0) cout << "got bbminmax_int: " << 
      //                     bbminmax_int[0] << ":" << -bbminmax_int[3] << " " << 
      //                     bbminmax_int[1] << ":" << -bbminmax_int[4] << " " << 
      //                     bbminmax_int[2] << ":" << -bbminmax_int[5] << endl;
      
      int my_count[7];
      int count[7];
      int split[3];
      int iter = 0;
      while (1) {

        // sorry this is so complicated. It is just the average of the 
        // bbox limits, but handles the integer rounding consistently 
        // for both negative and positive split's...
        FOR_I3 split[i] = bbminmax_int[i] - (bbminmax_int[i]+bbminmax_int[i+3])/2;
        
        FOR_I6 my_count[i] = 0;
        for (int ibb = 0; ibb < my_nbb; ++ibb) {
          // consider all 3 directions. We decide on the 
          // best one below...
          FOR_I3 {
            // note we use the smaller buf_double here for cache efficiency...
            if (buf_int[ibb*3+i] < split[i]) 
              ++my_count[i];
            else if (buf_int[ibb*3+i] == split[i])
              ++my_count[i+3];
          }
        }
        
        MPI_Allreduce(my_count,count,6,MPI_INT,MPI_SUM,mpi_split_comm);
        
        ++iter;
	if (!(iter < 100)) {
          cout << "XXXXXX " << mpi_rank << " mpi_split_rank: " << mpi_split_rank << endl;
        }
	assert(iter < 100);
      
        /*
        if (mpi_rank == 18) {
          //if (debug) {
          //cout << " bbminmax: 0: " << bbminmax[0] << ":" << -bbminmax[2] << 
          //" 1: " << bbminmax[1] << ":" << -bbminmax[3] << " split: " << split[0] << " " << split[1] << endl;
          cout << " DEBUG > " << iter << " got count[0,1,2]:  " << count[0] << " " << count[1] << " " << count[2] <<  
            " nbb/2: " << nbb_half << " nbb: " << nbb << " count[3,4,5]: " << count[3] << " " << count[4] << " " << count[5] << 
            " split[0,1,2]: " << split[0] << " " << split[1] << " " << split[2] << " bbminmax_int[2]: " << bbminmax_int[2] << " -bbminmax_int[2+3]: " << -bbminmax_int[2+3] << endl;
        }
        */
        
        // recall count[2] and count[3] contain the exact matches for directions j,k respectively...
        // this check determins if the matches are straddling the split...
        if ((count[0] <= nbb_half)&&(count[0]+count[3] >= nbb_half)&&
            (count[1] <= nbb_half)&&(count[1]+count[4] >= nbb_half)&&
            (count[2] <= nbb_half)&&(count[2]+count[5] >= nbb_half))
          break;
        
        // if we aren't done, do bisection...
        FOR_I3 {
          if (count[i]+count[i+3] < nbb_half) {
            bbminmax_int[i] = split[i];
          }
          else if (count[i] > nbb_half) {
            bbminmax_int[i+3] = -split[i];
          }
        }
          
      }
      
      // now everyone gets everyone's counts...
      my_count[6] = my_nbb;
      int (*my_count_of_rank)[7] = new int[mpi_split_size][7];
      MPI_Allgather(my_count,7,MPI_INT,(int*)my_count_of_rank,7,MPI_INT,mpi_split_comm);
      
      FOR_I7 count[i] = 0;
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        FOR_I3 count[i] += my_count_of_rank[rank][i];
        // recall [3,4,5] store the exact matches. These will be split between the 
        // 2 children so that the split is exactly: nbb/2, nbb-nbb/2
        if (rank < mpi_split_rank) {
          FOR_I3 count[3+i] += my_count_of_rank[rank][3+i]; // matches/equals for direction x,y,z
        }
        count[6] += my_count_of_rank[rank][6]; // totals -- same for both directions obviously
      }
      assert(count[6] == nbb); 
      
      const int mpi_split_size_half = mpi_split_size/2;
      
      // the left half will have nbb_half distributed as evenly as possible...
      int my_nbb_target0 = nbb_half/mpi_split_size_half;
      if (nbb_half%mpi_split_size_half) ++my_nbb_target0;
      int my_nbb_target1 = (nbb-nbb_half)/(mpi_split_size-mpi_split_size_half);
      if ((nbb-nbb_half)%(mpi_split_size-mpi_split_size_half)) ++my_nbb_target1;
      
      //if (mpi_rank == 0) cout << "GOT CHILD TARGETS: " << my_nbb_target0 << " " << my_nbb_target1 << endl;
      
      // We should be able to build the local/send pattern based on the above global reduction...
      // switch the my_count_of_rank[rank][0] to be the local count, and my_count_of_rank[rank][1] to be the send count
      
      int equal_count[3] = { 0, 0, 0 };
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        // my_count_of_rank[rank][0] contains the left count for rank. We may need to add some
        // or all of our equals...
        FOR_I3 {
          const int neq = max(0,min(my_count_of_rank[rank][3+i],nbb_half-count[i]-equal_count[i]));
          equal_count[i] += neq;
          my_count_of_rank[rank][i] += neq;
        }
      }
      
      double my_child_bbminmax[6][6];
      FOR_I6 FOR_J6 my_child_bbminmax[i][j] = HUGE_VAL; // max's are negative
      FOR_I6 my_count[i] = 0;
      FOR_I3 equal_count[i] = count[3+i]; // contains the count of equals up to but not including our rank
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        FOR_I3 {
          if (buf_int[ibb*3+i] < split[i]) {
            ++my_count[i];
            FOR_J6 my_child_bbminmax[2*i][j] = min(my_child_bbminmax[2*i][j],mybbVec[ibb].bbminmax[j]);
          }
          else if ((buf_int[ibb*3+i] == split[i])&&(count[i]+equal_count[i] < nbb_half)) {
            ++my_count[i];
            ++equal_count[i];
            FOR_J6 my_child_bbminmax[2*i][j] = min(my_child_bbminmax[2*i][j],mybbVec[ibb].bbminmax[j]);
          }
          else {
            FOR_J6 my_child_bbminmax[2*i+1][j] = min(my_child_bbminmax[2*i+1][j],mybbVec[ibb].bbminmax[j]);
          }
        }
      }
      FOR_I3 assert(my_count[i] == my_count_of_rank[mpi_split_rank][i]);
      
      double child_bbminmax[6][6];
      MPI_Allreduce((double*)my_child_bbminmax,(double*)child_bbminmax,36,MPI_DOUBLE,MPI_MIN,mpi_split_comm);
      double volume[3];
      FOR_I3 {
        const double volume0 = 
          (-child_bbminmax[2*i][3]-child_bbminmax[2*i][0])*
          (-child_bbminmax[2*i][4]-child_bbminmax[2*i][1])*
          (-child_bbminmax[2*i][5]-child_bbminmax[2*i][2]);
        assert(volume0 >= 0.0);
        const double volume1 = 
          (-child_bbminmax[2*i+1][3]-child_bbminmax[2*i+1][0])*
          (-child_bbminmax[2*i+1][4]-child_bbminmax[2*i+1][1])*
          (-child_bbminmax[2*i+1][5]-child_bbminmax[2*i+1][2]);
        assert(volume1 >= 0.0);
        volume[i] = volume0 + volume1;
      }
      
      // choose a direction for splitting based on minimizing the 
      // total child volume...
      
      int id = 0;
      if (volume[1] < volume[0]) {
        id = 1;
        if (volume[2] < volume[1])
          id = 2;
      }
      else if (volume[2] < volume[0]) {
        id = 2;
      }
      
      // =============================================
      // now with id selected, redistribute bboxes...
      // =============================================
      
      int * send_count = new int[mpi_split_size];
      int * recv_count = new int[mpi_split_size];
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        send_count[rank] = 0;
        recv_count[rank] = 0;
      }

      // do sends from 0 and recvs to 1 first...
      {
        int recv_rank = mpi_split_size_half-1; // start 1 less than the first 
        int nrecv = 0;
        for (int send_rank = 0; send_rank < mpi_split_size_half; ++send_rank) {
          // the number we need to send is...
          int nsend = my_count_of_rank[send_rank][6]-my_count_of_rank[send_rank][id];
          while (nsend > 0) {
            while (nrecv == 0) {
              // advance the recv_rank, and recompute nrecv
              ++recv_rank;
              assert(recv_rank < mpi_split_size);
              nrecv = max(0, my_nbb_target1 - (my_count_of_rank[recv_rank][6]-my_count_of_rank[recv_rank][id]) );
            }
            if (nrecv <= nsend) {
              // send_rank <-> recv_rank nrecv...
              if (send_rank == mpi_split_rank) send_count[recv_rank] += nrecv;
              if (recv_rank == mpi_split_rank) recv_count[send_rank] += nrecv;
              nsend -= nrecv;
              nrecv = 0;
            }
            else {
              // send_rank <-> recv_rank send_count...
              if (send_rank == mpi_split_rank) send_count[recv_rank] += nsend;
              if (recv_rank == mpi_split_rank) recv_count[send_rank] += nsend;
              nrecv -= nsend;
              nsend = 0;
            }
          }
        }
      }
        
      // then sends from 1 and recvs to 0...
      {
        int recv_rank = -1; // start one less than the first...
        int nrecv = 0;
        for (int send_rank = mpi_split_size_half; send_rank < mpi_split_size; ++send_rank) {
          // the number we need to send is...
          int nsend = my_count_of_rank[send_rank][id];
          while (nsend > 0) {
            while (nrecv == 0) {
              // advance the recv_rank, and recompute nrecv
              ++recv_rank;
              assert(recv_rank < mpi_split_size_half);
              nrecv = max(0, my_nbb_target0 - my_count_of_rank[recv_rank][id]);
            }
            if (nrecv <= nsend) {
              // send_rank <-> recv_rank nrecv...
              if (send_rank == mpi_split_rank) send_count[recv_rank] += nrecv;
              if (recv_rank == mpi_split_rank) recv_count[send_rank] += nrecv;
              nsend -= nrecv;
              nrecv = 0;
            }
            else {
              // send_rank <-> recv_rank send_count...
              if (send_rank == mpi_split_rank) send_count[recv_rank] += nsend;
              if (recv_rank == mpi_split_rank) recv_count[send_rank] += nsend;
              nrecv -= nsend;
              nsend = 0;
            }
          }
        }
      }

      delete[] my_count_of_rank;
        
      int * send_disp = new int[mpi_split_size];
      int * recv_disp = new int[mpi_split_size];

      /*
      // check of recv_count from send_count...
      MPI_Alltoall(send_count,1,MPI_INT,recv_disp,1,MPI_INT,mpi_split_comm);
      for (int rank = 0; rank < mpi_split_size; ++rank) 
      assert(recv_disp[rank] == recv_count[rank]);
      MPI_Barrier(mpi_split_comm);
      //if (mpi_rank == 0) cout << "recv_count looks good!" << endl;
      */

      send_disp[0] = 0;
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_split_size; ++rank) {
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
        recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];
      }
      const int send_count_sum = send_disp[mpi_split_size-1] + send_count[mpi_split_size-1];
      const int recv_count_sum = recv_disp[mpi_split_size-1] + recv_count[mpi_split_size-1];
        
      // compress the boxes that are staying local, and 
      // pack the ones that are going...

      //if (debug) cout << "DEBUG: DDDDDDDDDDDDDD id: " << id << " split[id]: " << split[id] << endl;
      //MPI_Pause("4");
      
      Bboxdata * send_buf = new Bboxdata[send_count_sum];
      int my_nbb_local = 0;
      int my_nbb_send = 0;
      int equal_count_id = count[3+id];
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        if (buf_int[ibb*3+id] < split[id]) {
          // this is a left bbox...
          if (mpi_split_rank < mpi_split_size_half) {
            // this bbox is where it belongs...
            mybbVec[my_nbb_local] = mybbVec[ibb];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 0 half...
            send_buf[my_nbb_send] = mybbVec[ibb];
            ++my_nbb_send;
          }
        }
        else if ((buf_int[ibb*3+id] == split[id])&&(count[id]+equal_count_id < nbb_half)) {
          ++equal_count_id;
          //if (debug) cout << "GOT AN EQUAL!!!!!!!!!!!!!!!!!" << endl;
          // this is also a left bbox...
          if (mpi_split_rank < mpi_split_size_half) {
            // this bbox is where it belongs...
            mybbVec[my_nbb_local] = mybbVec[ibb];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 0 half...
            send_buf[my_nbb_send] = mybbVec[ibb];
            ++my_nbb_send;
          }
        }
        else {
          // this is a right bbox...
          if (mpi_split_rank >= mpi_split_size_half) {
            // this bbox is where it belongs...
            mybbVec[my_nbb_local] = mybbVec[ibb];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 1 half...
            send_buf[my_nbb_send] = mybbVec[ibb];
            ++my_nbb_send;
          }
        }
      }
      assert(my_nbb_send == send_count_sum);

      delete[] buf_int;

      // resize the mybbVec arrays and do the exchange...
      int my_nbb_split = my_nbb_local + recv_count_sum;
      if (mybbVec.size() < my_nbb_split) mybbVec.resize(my_nbb_split);
      
      // check split counts...
      {
        int8 my_nbb_split_check[2];
        if (mpi_split_rank < mpi_split_size_half) {
          my_nbb_split_check[0] = my_nbb_split;
          my_nbb_split_check[1] = 0;
        }
        else {
          my_nbb_split_check[0] = 0;
          my_nbb_split_check[1] = my_nbb_split;
        }
        int8 nbb_split_check[2];
        MPI_Allreduce(my_nbb_split_check,nbb_split_check,2,MPI_INT8,MPI_SUM,mpi_split_comm);
        assert(nbb_split_check[0] == nbb_half);
        assert(nbb_split_check[1] == nbb-nbb_half);
      }

      // exchange bbox's directly into mybbVec. It was made large enough above...
      MPI_Alltoallv(send_buf,send_count,send_disp,MPI_Leafdata,
                    &mybbVec.front()+my_nbb_local,recv_count,recv_disp,MPI_Leafdata,
                    mpi_split_comm);
      
      delete[] send_buf;
      
      delete[] send_count;
      delete[] send_disp;
      delete[] recv_count;
      delete[] recv_disp;

      if (mpi_rank == 0) cout << " > time to bisect " << nbb << " on " << mpi_split_size << " ranks: " << MPI_Wtime()-wtime << endl;

      // now split and call this routine recursively...
        
      int mpi_key;
      if (mpi_split_rank < mpi_split_size_half)
        mpi_key = 0;
      else
        mpi_key = 1;
      
      MPI_Comm mpi_split_again_comm;
      MPI_Comm_split(mpi_split_comm, mpi_key, mpi_split_rank, &mpi_split_again_comm);
      
      // recursively call this routine...
      if (mpi_key == 0) {
        // this is the left...
        initRecursive(leafdata+leafdata_offset,mybbVec,my_nbb_split,nbb_half,child_bbminmax[2*id],2*ileaf+1,mpi_split_again_comm);
      }
      else {
        // this is the right...
        initRecursive(leafdata+leafdata_offset,mybbVec,my_nbb_split,nbb-nbb_half,child_bbminmax[2*id+1],2*ileaf+2,mpi_split_again_comm);
      }
      
      // and cleanup the split communicator...
      MPI_Comm_free(&mpi_split_again_comm);

    }
    else {

      assert(my_nbb != 0);
      assert(my_nbb == nbb);

      // figure out serial stack and leafdata sizes...
      int nstack = 1;
      int nleaves = 1;
      {
        int i = 1; 
        while (nleaves < my_nbb) {
          // nleaves needs the largest power of 2 greater than n...
          nleaves *= 2;
          // nstack needs 1 + 2 + 3 + 4 + 5 + etc...
          i += 1;
          nstack = nstack + i;
        }
        // and one final multiplier to allow all leaves to
        // contain just one member...
        nleaves *= 2;
        nleaves -= 1;
      }
      
      // we will need a stack.... 
      int * stack = new int[nstack];
      
      FOR_I2 assert(leafdata[0].index[i] == -2);
      
      // indicator leaf...
      leafdata[0].index[0] = -ileaf-1; // indicates a "special" leaf data...
      leafdata[0].index[1] = nleaves;

      // add the serial leaf data... 
      initSerial(leafdata+1,stack,mybbVec,nbb,bbminmax);
      
      // cleanup...
      delete[] stack;

    }
    
  }

  void initSerial(Leafdata *leafdata,int * stack,vector<Bboxdata>& mybbVec,const int my_nbb,double bbminmax[6]) {

    if (mpi_rank == 0) cout << " > starting serial build on my_nbb: " << my_nbb << endl;
    
    assert(my_nbb == nbb_serial_check);

    double wtime = MPI_Wtime();
    
    // check for a single bbox: i.e. terminal leaf...

    FOR_I2 assert(leafdata[0].index[i] == -2);

    if (my_nbb == 1) {
      
      FOR_I6 assert(bbminmax[i] == mybbVec[0].bbminmax[i]);
      FOR_I6 leafdata[0].bbminmax[i] = mybbVec[0].bbminmax[i];
      leafdata[0].index[0] = leafdata[0].index[1] = mybbVec[0].ist; // index
      return;

    }

    // otherwise set leaf 0...
    
    assert(my_nbb > 1);
    
    FOR_I6 leafdata[0].bbminmax[i] = bbminmax[i]; // jmin
    leafdata[0].index[0] = 0;
    leafdata[0].index[1] = my_nbb-1;

    stack[0] = 0;
    int istack = 1;

    // one sort...

    vector<pair<double,int> > idVec[3];
    FOR_I3 {
      idVec[i].resize(my_nbb);
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        idVec[i][ibb].first = mybbVec[ibb].bbminmax[i]; // xmin,ymin,zmin
        idVec[i][ibb].second = ibb;
      }
      sort(idVec[i].begin(),idVec[i].end());
    }

    if (mpi_rank == 0) cout << " > done sorts: " << MPI_Wtime()-wtime << endl;

    Bboxdata * bb_buf = new Bboxdata[my_nbb/2+1];
    pair<double,int> * pair_buf = new pair<double,int>[my_nbb/2+1]; // 5/2 ish
    int * flag = new int[my_nbb];
    while (istack > 0) {
      
      // pop the next leaf off the stack...
      const int ileaf = stack[--istack];
      //assert(2*ileaf+2 < nleaves);
      
      // the leaf data contains the range in mybbVec...
      const int ibb_left_f = leafdata[ileaf].index[0];
      const int ibb_right_l = leafdata[ileaf].index[1];
      
      // figure out the split...
      const int nbb = ibb_right_l-ibb_left_f+1;
      assert(nbb >= 2); // if it were 1 (never 0), should have been handled already!
      const int nbb_half = nbb/2;

      const int ibb_left_l = ibb_left_f + nbb_half - 1; 
      const int ibb_right_f = ibb_left_l + 1;

      // check child bounding boxes...
      double child_bbminmax[6][6];
      FOR_I6 FOR_J6 child_bbminmax[i][j] = HUGE_VAL;
      for (int id = 0; id <= 2; ++id) {
        for (int ibb_ = ibb_left_f; ibb_ <= ibb_left_l; ++ibb_) {
          const int ibb = idVec[id][ibb_].second;
          FOR_J6 child_bbminmax[2*id][j] = min(child_bbminmax[2*id][j],mybbVec[ibb].bbminmax[j]);
        }
        for (int ibb_ = ibb_right_f; ibb_ <= ibb_right_l; ++ibb_) {
          const int ibb = idVec[id][ibb_].second;
          FOR_J6 child_bbminmax[2*id+1][j] = min(child_bbminmax[2*id+1][j],mybbVec[ibb].bbminmax[j]);
        }
      }
      
      double volume[3];
      FOR_I3 {
        const double volume0 = 
          (-child_bbminmax[2*i][3]-child_bbminmax[2*i][0])*
          (-child_bbminmax[2*i][4]-child_bbminmax[2*i][1])*
          (-child_bbminmax[2*i][5]-child_bbminmax[2*i][2]);
        assert(volume0 >= 0.0);
        const double volume1 = 
          (-child_bbminmax[2*i+1][3]-child_bbminmax[2*i+1][0])*
          (-child_bbminmax[2*i+1][4]-child_bbminmax[2*i+1][1])*
          (-child_bbminmax[2*i+1][5]-child_bbminmax[2*i+1][2]);
        assert(volume1 >= 0.0);
        volume[i] = volume0 + volume1;
      }
      
      // choose a direction for splitting based on minimizing the 
      // total child volume...
      
      int id = 0;
      if (volume[1] < volume[0]) {
        id = 1;
        if (volume[2] < volume[1])
          id = 2;
      }
      else if (volume[2] < volume[0]) {
        id = 2;
      }

      //set the flag to 0, or 1, depending on which group we are going to...
      
      for (int ibb_ = ibb_left_f; ibb_ <= ibb_left_l; ++ibb_) {
        const int ibb = idVec[id][ibb_].second;
        flag[ibb] = 0;
      }
      for (int ibb_ = ibb_right_f; ibb_ <= ibb_right_l; ++ibb_) {
        const int ibb = idVec[id][ibb_].second;
        flag[ibb] = 1;
      }

      // reorder mybbVec[ibb_left_f:ibb_right_l] to respect the new division...
      int nbb_left = 0;
      int nbb_right = 0;
      for (int ibb = ibb_left_f; ibb <= ibb_right_l; ++ibb) {
        if (flag[ibb] == 0) {
          mybbVec[ibb_left_f+nbb_left] = mybbVec[ibb];
          flag[ibb] = ibb_left_f+nbb_left;
          ++nbb_left;
        }
        else {
          assert(flag[ibb] == 1);
          assert(nbb_right < my_nbb/2+1);
          bb_buf[nbb_right] = mybbVec[ibb];
          flag[ibb] = ibb_right_f+nbb_right;
          ++nbb_right;
        }
      }
      assert(nbb_left == nbb_half); 
      assert(nbb_right == nbb-nbb_half);

      for (int ibb = 0; ibb < nbb_right; ++ibb)
        mybbVec[ibb_right_f+ibb] = bb_buf[ibb];
      
      // for the idVec[id], it is already in the correct order, and we can just change its ibb (i.e. second)
      // reference...
      for (int ibb_ = ibb_left_f; ibb_ <= ibb_right_l; ++ibb_) {
        const int ibb = idVec[id][ibb_].second;
        idVec[id][ibb_].second = flag[ibb];
      }
      
      // and the other directions 1,2...
      for (int i = 1; i <= 2; ++i) {
        const int id_other = (id+i)%3;
        nbb_left = 0;
        nbb_right = 0;
        for (int ibb_ = ibb_left_f; ibb_ <= ibb_right_l; ++ibb_) {
          const int ibb_old = idVec[id_other][ibb_].second;
          const int ibb = flag[ibb_old];
          if (ibb <= ibb_left_l) {
            // we are in the left part...
            idVec[id_other][ibb_left_f+nbb_left].first = idVec[id_other][ibb_].first;
            idVec[id_other][ibb_left_f+nbb_left].second = ibb;
            ++nbb_left;
          }
          else {
            assert(nbb_right < my_nbb/2+1); 
            pair_buf[nbb_right].first = idVec[id_other][ibb_].first;
            pair_buf[nbb_right].second = ibb;
            ++nbb_right;
          }
        }
        for (int ibb = 0; ibb < nbb_right; ++ibb) {
          idVec[id_other][ibb_right_f+ibb] = pair_buf[ibb];
        }
      }

      // now push these new leaves onto the stack...
      
      // left...

      FOR_I2 assert(leafdata[2*ileaf+1].index[i] == -2);
      FOR_I6 leafdata[2*ileaf+1].bbminmax[i] = child_bbminmax[2*id][i]; 
      
      if (nbb_left == 1) {
        // this is the terminal leaf...
        leafdata[2*ileaf+1].index[0] = leafdata[2*ileaf+1].index[1] = mybbVec[ibb_left_f].ist; // index
      }
      else {
        // put the leaf on the stack...
        assert(nbb_left > 1);
        leafdata[2*ileaf+1].index[0] = ibb_left_f;
        leafdata[2*ileaf+1].index[1] = ibb_left_l;
        stack[istack++] = 2*ileaf+1;
      }
      
      // right...
      
      FOR_I2 assert(leafdata[2*ileaf+2].index[i] == -2);
      FOR_I6 leafdata[2*ileaf+2].bbminmax[i] = child_bbminmax[2*id+1][i];
      
      if (nbb_right == 1) {
        // this is the terminal leaf...
        leafdata[2*ileaf+2].index[0] = leafdata[2*ileaf+2].index[1] = mybbVec[ibb_right_l].ist; // index
      }
      else {
        // put the leaf on the stack...
        assert(nbb_right > 1);
        leafdata[2*ileaf+2].index[0] = ibb_right_f;
        leafdata[2*ileaf+2].index[1] = ibb_right_l;
        stack[istack++] = 2*ileaf+2;
      }
      
    }
    
    delete[] flag;
    delete[] pair_buf;
    delete[] bb_buf;
    
    if (mpi_rank == 0) cout << " > done serial: " << MPI_Wtime()-wtime << endl;
    
  }

};

#endif















#ifdef JUNKJUNK

    // the "new" method fails on large problems, probably for the same reason we
    // had to split up the bcast below...

    // =============================================================
    // gather to rank0 of shm, then gatherallv...
    // =============================================================
    
    int * nleaves_of_rank = NULL;
    if (mpi_rank_shared == 0) nleaves_of_rank = new int[mpi_size_shared];
    MPI_Gather(&my_nleaves,1,MPI_INT,nleaves_of_rank,1,MPI_INT,0,mpi_comm_shared);

    int * disp_of_rank = NULL;
    Leafdata * leafdata_buf_shared = NULL;
    int leafdata_count_shared;
    if (mpi_rank_shared == 0) {
      // check for overflow...
      int8 leafdata_count_int8 = 0;
      for (int rank = 0; rank < mpi_size_shared; ++rank) leafdata_count_int8 += nleaves_of_rank[rank];
      assert(leafdata_count_int8 < TWO_BILLION);
      disp_of_rank = new int[mpi_size_shared];
      disp_of_rank[0] = 0;
      for (int rank = 1; rank < mpi_size_shared; ++rank)
        disp_of_rank[rank] = disp_of_rank[rank-1] + nleaves_of_rank[rank-1];
      leafdata_count_shared = disp_of_rank[mpi_size_shared-1] + nleaves_of_rank[mpi_size_shared-1];
      leafdata_buf_shared = new Leafdata[leafdata_count_shared];
    }
    
    MPI_Gatherv(my_leafdata,my_nleaves,MPI_Leafdata,leafdata_buf_shared,nleaves_of_rank,disp_of_rank,MPI_Leafdata,0,mpi_comm_shared);
    delete[] my_leafdata;
    
    if (mpi_rank == 0) {
      delete[] nleaves_of_rank;
      delete[] disp_of_rank;
    }
    
    if (mpi_rank == 0) cout << " > done gather to all shared: " << MPI_Wtime()-wtime0 << endl;
    
    // at this point, all the ranks of a given shared node have there leaf data
    // gathered at mpi_rank_shared == 0 in leafdata_buf_shared.

    // now gatherallv to all rank 0's using the internode comm...
    
    Leafdata * leafdata_buf = NULL;
    int leafdata_count;
    
    if (mpi_rank_shared == 0) {
      
      nleaves_of_rank = new int[mpi_size_internode];
      MPI_Allgather(&leafdata_count_shared,1,MPI_INT,nleaves_of_rank,1,MPI_INT,mpi_comm_internode);

      int8 leafdata_count_int8 = 0;
      for (int rank = 0; rank < mpi_size_internode; ++rank) leafdata_count_int8 += nleaves_of_rank[rank];
      assert(leafdata_count_int8 < TWO_BILLION);
      disp_of_rank = new int[mpi_size_internode];
      disp_of_rank[0] = 0;
      for (int rank = 1; rank < mpi_size_internode; ++rank)
        disp_of_rank[rank] = disp_of_rank[rank-1] + nleaves_of_rank[rank-1];
      leafdata_count = disp_of_rank[mpi_size_internode-1] + nleaves_of_rank[mpi_size_internode-1];
      leafdata_buf = new Leafdata[leafdata_count];
      
      MPI_Allgatherv(leafdata_buf_shared,leafdata_count_shared,MPI_Leafdata,leafdata_buf,nleaves_of_rank,disp_of_rank,MPI_Leafdata,mpi_comm_internode);
    
      delete[] leafdata_buf_shared;
      delete[] nleaves_of_rank;
      delete[] disp_of_rank;
    }

    // everyone allocates the leafdata in shared memory...
    assert(leafdata == NULL);
    CTI_Mmap(leafdata,nleaves);
    
    if (mpi_rank == 0) cout << " > done shm allocation: " << MPI_Wtime()-wtime0 << endl;
    
    // to avoid the bcast, all mpi_rank_shared == 0 unpacks into the shared memory...
    
    if (mpi_rank_shared == 0) {
      
      // for checking...
      //for (int ileaf = 0; ileaf < nleaves; ++ileaf)
      //  leafdata[ileaf].index[0] = -1;

      // unpack the leafdata_buf into the shared memory leafdata...
      int ileafdata = 0;
      while (ileafdata < leafdata_count) {
        
        // the first leaf from any rank should be a "special" leaf 
        // that tells us how to unpack (use -1 indexing)... 
        assert(leafdata_buf[ileafdata].index[0] < 0);
        int ileaf = -leafdata_buf[ileafdata].index[0]-1;
        int nleaf = leafdata_buf[ileafdata].index[1]; // nleaf = 1,3,7,15,...2^x-1
        ++ileafdata;
        
        int nlevel = 1; // 1,2,4,8,16,...
        while (nleaf > 0) {
          // pack nlevel leafdata elements into the ileaf location... 
          for (int ilevel = 0; ilevel < nlevel; ++ilevel) {
            assert(ileaf+ilevel < nleaves);
            //assert(leafdata[ileaf+ilevel].index[0] == -1);
            leafdata[ileaf+ilevel] = leafdata_buf[ileafdata+ilevel];
          }
          // get ready for the next level...
	  ileafdata += nlevel;
          nleaf -= nlevel;
          nlevel *= 2;
          ileaf = 2*ileaf+1;
        }
        // we should have gotten to exactly 0...
        assert(nleaf == 0);
        
      }
      
      delete[] leafdata_buf;
      
    }
    
    if (mpi_rank == 0) cout << " > done repacking: " << MPI_Wtime()-wtime0 << endl;

#endif
