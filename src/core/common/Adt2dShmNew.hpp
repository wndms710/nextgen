#ifndef _ADT2DSHMNEW_HPP_
#define _ADT2DSHMNEW_HPP_

#include <stack>

class Adt2dShmNew {
private:

  typedef struct {
    int bbminmax[4];
    int index[2];
  } Leafdata;
  
  int nbb_serial_check;
  double wtime0;
  int nleaves;
  Leafdata * leafdata;
  
#define T Leafdata
# include "Mmap_Munmap.hpp"
#undef T

public:

  Adt2dShmNew(vector<int>& mybbVec,const int my_nbb) {
    
    if (mpi_rank == 0) cout << "Adt2dShmNew()" << endl;
    
    wtime0 = MPI_Wtime();
    nbb_serial_check = -1;
    nleaves = 0;
    leafdata = NULL;
    
    assert(mybbVec.size() >= 5*my_nbb);

    // define the Leafdata struct...
    assert(int_size   == 4ll);
    assert(sizeof(Leafdata) == 24);
    
    MPI_Datatype oldtypes[2];
    int blockcounts[2];
    MPI_Aint offsets[2];
    
    // recall Leafdata...
    //typedef struct {
    //  int bbminmax[4];
    //  int index[2];
    //} Leafdata;
    
    // int bbminmax[4]...
    offsets[0] = 0;
    oldtypes[0] = MPI_INT;
    blockcounts[0] = 4;

    // int index[2]...
    offsets[1] = offsets[0] + int_size*4;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 2;
    
    // build the Leafdata type...
    MPI_Datatype MPI_Leafdata;
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &MPI_Leafdata);
    MPI_Type_commit(&MPI_Leafdata);

    // check that the type size matches the struct size...
    int leafdata_size;
    MPI_Type_size(MPI_Leafdata,&leafdata_size);
    assert(leafdata_size == sizeof(Leafdata));
    
    // Compute the global bounding box...
    int my_bbminmax[4] = { TWO_BILLION, TWO_BILLION, TWO_BILLION, TWO_BILLION };
    for (int ibb = 0; ibb < my_nbb; ++ibb) {
      FOR_I4 my_bbminmax[i] = min(my_bbminmax[i],mybbVec[5*ibb+1+i]);
    }
    int bbminmax[4];
    MPI_Allreduce(my_bbminmax,bbminmax,4,MPI_INT,MPI_MIN,mpi_comm);
    
    // and the global count...
    int8 my_nbb_int8 = my_nbb;
    int8 nbb_int8;
    MPI_Allreduce(&my_nbb_int8,&nbb_int8,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert(nbb_int8 < TWO_BILLION);
    const int nbb = nbb_int8;
    
    if (mpi_rank == 0) cout << " > got global nbb: " << nbb << " with bbox: " << bbminmax[0] << " " << bbminmax[1] << " " << bbminmax[2] << " " << bbminmax[3] << endl;
    
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
    
    // figure out how big our part of the leaf data will be...
    int my_nleaves = 0;
    int my_nbb_split = nbb;
    {
      int my_size_split = mpi_size;
      int my_rank_split = mpi_rank;
      while (my_size_split > 1) {
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
    int (*my_leafdata)[6] = new int[my_nleaves][6];
    
    // for now, fill with -2 for checking...
    for (int ileaf = 0; ileaf < my_nleaves; ++ileaf)
      FOR_I6 my_leafdata[ileaf][i] = -2;
    
    if (mpi_rank == 0) cout << " > about to call initRecursive..." << endl;
    
    // now call the recursive split and shuffle alorithm...
    initRecursive(my_leafdata,mybbVec,my_nbb,nbb,bbminmax,0,mpi_comm);

    // now gather the leaf data at global rank 0...
    int * nleaves_of_rank = NULL;
    if (mpi_rank == 0) nleaves_of_rank = new int[mpi_size];
    MPI_Gather(&my_nleaves,1,MPI_INT,nleaves_of_rank,1,MPI_INT,0,mpi_comm);
    
    int * disp_of_rank = NULL;
    Leafdata * leafdata_buf = NULL;
    int leafdata_count;
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
    
    MPI_Gatherv((Leafdata*)my_leafdata,my_nleaves,MPI_Leafdata,leafdata_buf,nleaves_of_rank,disp_of_rank,MPI_Leafdata,0,mpi_comm);
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
        // that tells us how to unpack... 
        assert(leafdata_buf[ileafdata].index[0] == -1);
        int ileaf = leafdata_buf[ileafdata].bbminmax[0];
        int nleaf = leafdata_buf[ileafdata].bbminmax[1]; // nleaf = 1,3,7,15,...2^x-1
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

  ~Adt2dShmNew() {
    
    // free the shm...
    CTI_Munmap(leafdata,nleaves);
    
  }

  // public check function...
  
  void getRandomIstAndBbox(int& ist,int bbmin[2],int bbmax[2]) {
    
    int ileaf = 0;
    // traverse to the bottom...
    while (leafdata[ileaf].index[0] != leafdata[ileaf].index[1]) {
      // choose a random direction...
      int id = rand()%2; // 0 or 1
      ileaf = 2*ileaf+1+id;
    }
    ist = leafdata[ileaf].index[0];
    FOR_I2 bbmin[i] = leafdata[ileaf].bbminmax[i];
    FOR_I2 bbmax[i] = leafdata[ileaf].bbminmax[i+2];
    
  }

  void getGlobalBbox(int bbmin[2],int bbmax[2]) const {

    // the global bbox is the bbox of the first leaf...
    assert(nleaves > 0);
    FOR_I2 bbmin[i] = leafdata[0].bbminmax[i];
    FOR_I2 bbmax[i] = leafdata[0].bbminmax[i+2];
    
  }
  
  void addBboxesForPoint(vector<int>& bboxVec,const int jk[2]) const {
    
    if (nleaves > 0) {

      stack<int> stack;
      stack.push(0);

      while (!stack.empty()) {
      
        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();

        // skip if out of bounds...
        if ( (jk[0] > leafdata[ileaf].bbminmax[2]) || (jk[0] < leafdata[ileaf].bbminmax[0]) ||
             (jk[1] > leafdata[ileaf].bbminmax[3]) || (jk[1] < leafdata[ileaf].bbminmax[1]) )
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

  void addBboxesForCircle(vector<int>& bboxVec,const int jk[2],const int rjk) const {
    
    if (nleaves > 0) {

      stack<int> stack;
      stack.push(0);
      
      const int8 r2jk = int8(rjk)*int8(rjk);
      while (!stack.empty()) {
      
        // pop the next leaf off the stack...
        const int ileaf = stack.top(); stack.pop();

	// calculate the distance of jk to the bbox...
	int8 d2 = 0;
	FOR_I2 {
	  int8 djk = max(0,max(leafdata[ileaf].bbminmax[i]-jk[i],jk[i]-leafdata[ileaf].bbminmax[2+i]));
	  d2 += djk*djk;
	}

        // skip if out of bounds...
        if (d2 >= r2jk) 
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
   
private:

  void initRecursive(int (*leafdata)[6],vector<int>& mybbVec,const int my_nbb,const int nbb,int bbminmax[4],const int ileaf,const MPI_Comm& mpi_split_comm) {
    
    int mpi_split_size;
    MPI_Comm_size(mpi_split_comm,&mpi_split_size);
      
    if (mpi_split_size > 1) {

      double wtime = MPI_Wtime();
      
      // reduce to get nbb...
      const int nbb_half = nbb/2;

      int mpi_split_rank;
      MPI_Comm_rank(mpi_split_comm, &mpi_split_rank);
      
      // if we are rank0 on this set of data, we are responsible for adding its
      // info to the leaf data we will be returning to master...

      int leafdata_offset = 0;
      if (mpi_split_rank == 0) {
        
        // used below...
        leafdata_offset = 2;
        
        FOR_I6 assert(leafdata[0][i] == -2);
        
        // indicator leaf...
        leafdata[0][0] = ileaf;
        leafdata[0][1] = 1;
        leafdata[0][4] = -1; // indicates a "special" leaf data...

        FOR_I6 assert(leafdata[1][i] == -2);
        
        // bbox leaf...
        leafdata[1][0] = bbminmax[0];
        leafdata[1][1] = bbminmax[1];
        leafdata[1][2] = -bbminmax[2];
        leafdata[1][3] = -bbminmax[3];
        leafdata[1][4] = 0; 
        leafdata[1][5] = my_nbb-1; 
        
      }
      
      // compute the global bounding box, and also store the minimum bounds of all boxes
      // for cache efficiency...
        
      int * buf_int = new int[my_nbb*2];
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        FOR_I2 buf_int[2*ibb+i] = mybbVec[5*ibb+1+i]; // j min, k min
      }
        
      // now would be a good time to store this bbminmax in the parent leaf. It is 
      // about to be modified...

      int my_count[5];
      int count[5];
      int split[2];
      int iter = 0;
      while (1) {

        // sorry this is so compicated. It is just the average of the 
        // bbox limits, but handles the integer rounding consistently 
        // for both negative and positive split's...
        FOR_I2 split[i] = bbminmax[i] - (bbminmax[i]+bbminmax[i+2])/2;

        FOR_I4 my_count[i] = 0;
        for (int ibb = 0; ibb < my_nbb; ++ibb) {
          // consider both directions. We decide on the 
          // best one below...
          FOR_I2 {
            // note we use the smaller buf_int here for cache efficiency...
            if (buf_int[2*ibb+i] < split[i]) 
              ++my_count[i];
            else if (buf_int[2*ibb+i] == split[i])
              ++my_count[i+2];
          }
        }
        
        MPI_Allreduce(my_count,count,4,MPI_INT,MPI_SUM,mpi_split_comm);
        
        ++iter;
	assert(iter < 100);
        /*
          if (mpi_rank == 0) {
          //if (debug) {
          //cout << " bbminmax: 0: " << bbminmax[0] << ":" << -bbminmax[2] << 
          //" 1: " << bbminmax[1] << ":" << -bbminmax[3] << " split: " << split[0] << " " << split[1] << endl;
          cout << " DEBUG > " << iter << " got count[0]:  " << count[0] << " count[1]: " << count[1] << 
          " nbb/2: " << nbb_half << " nbb: " << nbb << " count[2]: " << count[2] << " count[3]: " << count[3] << endl;
          }
        */

        // recall count[2] and count[3] contain the exact matches for directions j,k respectively...
        // this check determins if the matches are straddling the split...
        if ((count[0] <= nbb_half)&&(count[0]+count[2] >= nbb_half)&&
            (count[1] <= nbb_half)&&(count[1]+count[3] >= nbb_half))
          break;
        
        // if we aren't done, do bisection...
        /*
          FOR_I2 {
          if (count[i] < nbb_half) {
          bbminmax[i] = split[i];
          }
          else {
          bbminmax[i+2] = -split[i];
          }
          }
        */
        // I think this is better...
        FOR_I2 {
          if (count[i]+count[i+2] < nbb_half) {
            bbminmax[i] = split[i];
          }
          else if (count[i] > nbb_half) {
            bbminmax[i+2] = -split[i];
          }
        }
          
      }
      
      delete[] buf_int;

      // now everyone gets everyone's counts...
      my_count[4] = my_nbb;
      int (*my_count_of_rank)[5] = new int[mpi_split_size][5];
      MPI_Allgather(my_count,5,MPI_INT,(int*)my_count_of_rank,5,MPI_INT,mpi_split_comm);
      
      FOR_I5 count[i] = 0;
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        FOR_I2 count[i] += my_count_of_rank[rank][i];
        // recall [2],[3] store the exact matches. These will be split between the 
        // 2 children so that the split is exactly: nbb/2, nbb-nbb/2
        if (rank < mpi_split_rank) {
          FOR_I2 count[2+i] += my_count_of_rank[rank][2+i]; // matches/equals for direction j,k
        }
        count[4] += my_count_of_rank[rank][4]; // totals -- same for both directions obviously
      }
      assert(count[4] == nbb); 

      const int mpi_split_size_half = mpi_split_size/2;

      // the left half will have nbb_half distributed as evenly as possible...
      int my_nbb_target0 = nbb_half/mpi_split_size_half;
      if (nbb_half%mpi_split_size_half) ++my_nbb_target0;
      int my_nbb_target1 = (nbb-nbb_half)/(mpi_split_size-mpi_split_size_half);
      if ((nbb-nbb_half)%(mpi_split_size-mpi_split_size_half)) ++my_nbb_target1;

      //if (mpi_rank == 0) cout << "GOT CHILD TARGETS: " << my_nbb_target0 << " " << my_nbb_target1 << endl;

      // We should be able to build the local/send pattern based on the above global reduction...
      // switch the my_count_of_rank[rank][0] to be the local count, and my_count_of_rank[rank][1] to be the send count

      int equal_count[2] = { 0, 0 };
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        // my_count_of_rank[rank][0] contains the left count for rank. We may need to add some
        // or all of our equals...
        FOR_I2 {
          const int neq = max(0,min(my_count_of_rank[rank][2+i],nbb_half-count[i]-equal_count[i]));
          equal_count[i] += neq;
          my_count_of_rank[rank][i] += neq;
        }
      }
        
      int my_child_bbminmax[4][4];
      FOR_I4 FOR_J4 my_child_bbminmax[i][j] = TWO_BILLION; // max's are negative
      FOR_I4 my_count[i] = 0;
      FOR_I2 equal_count[i] = count[2+i]; // contains the count of equals up to but not including our rank
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        FOR_I2 {
          if (mybbVec[5*ibb+1+i] < split[i]) {
            ++my_count[i];
            FOR_J4 my_child_bbminmax[2*i][j] = min(my_child_bbminmax[2*i][j],mybbVec[5*ibb+1+j]);
          }
          else if ((mybbVec[5*ibb+1+i] == split[i])&&(count[i]+equal_count[i] < nbb_half)) {
            ++my_count[i];
            ++equal_count[i];
            FOR_J4 my_child_bbminmax[2*i][j] = min(my_child_bbminmax[2*i][j],mybbVec[5*ibb+1+j]);
          }
          else {
            FOR_J4 my_child_bbminmax[2*i+1][j] = min(my_child_bbminmax[2*i+1][j],mybbVec[5*ibb+1+j]);
          }
        }
      }
      FOR_I2 assert(my_count[i] == my_count_of_rank[mpi_split_rank][i]);
      
      int child_bbminmax[4][4];
      MPI_Allreduce((int*)my_child_bbminmax,(int*)child_bbminmax,16,MPI_INT,MPI_MIN,mpi_split_comm);
      int8 area[2];
      FOR_I2 {
        const int8 area0 = 
          int8(child_bbminmax[2*i][0]  +child_bbminmax[2*i][2])*
          int8(child_bbminmax[2*i][1]  +child_bbminmax[2*i][3]);
        const int8 area1 = 
          int8(child_bbminmax[2*i+1][0]+child_bbminmax[2*i+1][2])*
          int8(child_bbminmax[2*i+1][1]+child_bbminmax[2*i+1][3]);
        area[i] = area0 + area1;
      }

      // choose a direction for splitting based on minimizing the 
      // child area...
    
      int id = 0;
      if (area[1] < area[0])
        id = 1;
      
      //if (debug) cout << "DEBUG: DDDDDDDDDDDDDD id: " << id << " split[id]: " << split[id] << endl;
      
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
          int nsend = my_count_of_rank[send_rank][4]-my_count_of_rank[send_rank][id];
          while (nsend > 0) {
            while (nrecv == 0) {
              // advance the recv_rank, and recompute nrecv
              ++recv_rank;
              assert(recv_rank < mpi_split_size);
              nrecv = max(0, my_nbb_target1 - (my_count_of_rank[recv_rank][4]-my_count_of_rank[recv_rank][id]) );
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
      
      int * send_buf_int = new int[5*send_count_sum];
      int my_nbb_local = 0;
      int my_nbb_send = 0;
      int equal_count_id = count[2+id];
      for (int ibb = 0; ibb < my_nbb; ++ibb) {
        if (mybbVec[5*ibb+1+id] < split[id]) {
          // this is a left bbox...
          if (mpi_split_rank < mpi_split_size_half) {
            // this bbox is where it belongs...
            FOR_I5 mybbVec[5*my_nbb_local+i] = mybbVec[5*ibb+i];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 0 half...
            FOR_I5 send_buf_int[5*my_nbb_send+i] = mybbVec[5*ibb+i];
            ++my_nbb_send;
          }
        }
        else if ((mybbVec[5*ibb+1+id] == split[id])&&(count[id]+equal_count_id < nbb_half)) {
          ++equal_count_id;
          //if (debug) cout << "GOT AN EQUAL!!!!!!!!!!!!!!!!!" << endl;
          // this is also a left bbox...
          if (mpi_split_rank < mpi_split_size_half) {
            // this bbox is where it belongs...
            FOR_I5 mybbVec[5*my_nbb_local+i] = mybbVec[5*ibb+i];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 0 half...
            FOR_I5 send_buf_int[5*my_nbb_send+i] = mybbVec[5*ibb+i];
            ++my_nbb_send;
          }
        }
        else {
          // this is a right bbox...
          if (mpi_split_rank >= mpi_split_size_half) {
            // this bbox is where it belongs...
            FOR_I5 mybbVec[5*my_nbb_local+i] = mybbVec[5*ibb+i];
            ++my_nbb_local;
          }
          else {
            // this bbox must be sent to the 1 half...
            FOR_I5 send_buf_int[5*my_nbb_send+i] = mybbVec[5*ibb+i];
            ++my_nbb_send;
          }
        }
      }
      assert(my_nbb_send == send_count_sum);

      // resize the mybbVec arrays and do the exchange...
      int my_nbb_split = my_nbb_local + recv_count_sum;
      if (mybbVec.size() < 5*my_nbb_split) mybbVec.resize(5*my_nbb_split);
      
      /*
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
      */
        
      for (int rank = 0; rank < mpi_split_size; ++rank) {
        send_count[rank] *= 5;
        recv_count[rank] *= 5;
        send_disp[rank] *= 5;
        recv_disp[rank] *= 5;
      }

      // exchange bbox's directly into mybbVec. It was made large enough above...
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                    &mybbVec.front()+my_nbb_local*5,recv_count,recv_disp,MPI_INT,
                    mpi_split_comm);
      
      delete[] send_buf_int;
      
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

      // 
      int * stack = new int[nstack];
      
      FOR_I6 assert(leafdata[0][i] == -2);
      
      // indicator leaf...
      leafdata[0][0] = ileaf;
      leafdata[0][1] = nleaves;
      leafdata[0][4] = -1; // indicates a "special" leaf data...

      //  
      initSerial(leafdata+1,stack,mybbVec,nbb,bbminmax);
      
      delete[] stack;
      
    }
    
  }
  
  void initSerial(int (*leafdata)[6],int * stack,vector<int>& mybbVec,const int nbb,int bbminmax[4]) {

    if (mpi_rank == 0) cout << " > starting serial build..." << endl;
    
    assert(nbb == nbb_serial_check);

    double wtime = MPI_Wtime();
    
    // set leaf 0...
    
    // check...
    FOR_I6 assert(leafdata[0][i] == -2);

    leafdata[0][0] = bbminmax[0]; // jmin
    leafdata[0][1] = bbminmax[1]; // kmin
    leafdata[0][2] = -bbminmax[2]; assert(leafdata[0][2] >= leafdata[0][0]);
    leafdata[0][3] = -bbminmax[3]; assert(leafdata[0][3] >= leafdata[0][1]);
    leafdata[0][4] = 0;
    leafdata[0][5] = nbb-1;
    
    if (nbb == 1) {

      // this is the terminal leaf...
      leafdata[0][4] = leafdata[0][5] = mybbVec[0]; // index
      
    }
    else {
      
      assert(nbb > 1);
    
      stack[0] = 0;
      int istack = 1;

      vector<pair<int,int> > idVec[2];

      idVec[0].resize(nbb);
      for (int ibb = 0; ibb < nbb; ++ibb) {
        idVec[0][ibb].first = mybbVec[5*ibb+1]; // jmin
        idVec[0][ibb].second = ibb;
      }
      sort(idVec[0].begin(),idVec[0].end());
    
      idVec[1].resize(nbb);
      for (int ibb = 0; ibb < nbb; ++ibb) {
        idVec[1][ibb].first = mybbVec[5*ibb+2]; // kmin
        idVec[1][ibb].second = ibb;
      }
      sort(idVec[1].begin(),idVec[1].end());

      if (mpi_rank == 0) cout << " > done sorts: " << MPI_Wtime()-wtime << endl;
    
      int * int_buf = new int[nbb*3+1]; // 5/2 ish
      int * flag = new int[nbb];
      while (istack > 0) {
      
        // pop the next leaf off the stack...
        const int ileaf = stack[--istack];
        //assert(2*ileaf+2 < nleaves);
        
        // the leaf data contains the range in mybbVec...
        const int ibb_left_f = leafdata[ileaf][4];
        const int ibb_right_l = leafdata[ileaf][5];
      
        // figure out the split...
        const int nbb = ibb_right_l-ibb_left_f+1;
        assert(nbb >= 2); // if it were 1 (never 0), should have been handled already!
        const int nbb_half = nbb/2;

        const int ibb_left_l = ibb_left_f + nbb_half - 1; 
        const int ibb_right_f = ibb_left_l + 1;

        // check child bounding boxes...
        int child_bbminmax[4][4];
        FOR_I4 FOR_J4 child_bbminmax[i][j] = TWO_BILLION; // max's are negative
        for (int id = 0; id <= 1; ++id) {
          for (int ibb_ = ibb_left_f; ibb_ <= ibb_left_l; ++ibb_) {
            const int ibb = idVec[id][ibb_].second;
            FOR_J4 child_bbminmax[2*id][j] = min(child_bbminmax[2*id][j],mybbVec[5*ibb+1+j]);
          }
          for (int ibb_ = ibb_right_f; ibb_ <= ibb_right_l; ++ibb_) {
            const int ibb = idVec[id][ibb_].second;
            FOR_J4 child_bbminmax[2*id+1][j] = min(child_bbminmax[2*id+1][j],mybbVec[5*ibb+1+j]);
          }
        }
      
        int8 area[2];
        FOR_I2 {
          const int8 area0 = 
            int8(child_bbminmax[2*i][0]  +child_bbminmax[2*i][2])*
            int8(child_bbminmax[2*i][1]  +child_bbminmax[2*i][3]);
          const int8 area1 = 
            int8(child_bbminmax[2*i+1][0]+child_bbminmax[2*i+1][2])*
            int8(child_bbminmax[2*i+1][1]+child_bbminmax[2*i+1][3]);
          area[i] = area0 + area1;
        }
      
        // choose a direction for splitting based on minimizing the 
        // child area...
        int id = 0;
        if (area[1] < area[0])
          id = 1;
      
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
            FOR_I5 mybbVec[5*(ibb_left_f+nbb_left)+i] = mybbVec[5*ibb+i];
            flag[ibb] = ibb_left_f+nbb_left;
            ++nbb_left;
          }
          else {
            FOR_I5 int_buf[5*nbb_right+i] = mybbVec[5*ibb+i];
            flag[ibb] = ibb_right_f+nbb_right;
            ++nbb_right;
          }
        }
        assert(nbb_left == nbb_half); 
        assert(nbb_right == nbb-nbb_half);

        for (int ibb = 0; ibb < nbb_right; ++ibb)
          FOR_I5 mybbVec[5*(ibb_right_f+ibb)+i] = int_buf[5*ibb+i];
      
        // for the idVec[id], it is already in the correct order, and we can just change its ibb (i.e. second)
        // reference...
        for (int ibb_ = ibb_left_f; ibb_ <= ibb_right_l; ++ibb_) {
          const int ibb = idVec[id][ibb_].second;
          idVec[id][ibb_].second = flag[ibb];
        }
      
        // and the other direction...
        nbb_left = 0;
        nbb_right = 0;
        for (int ibb_ = ibb_left_f; ibb_ <= ibb_right_l; ++ibb_) {
          const int ibb_old = idVec[1-id][ibb_].second;
          const int ibb = flag[ibb_old];
          if (ibb <= ibb_left_l) {
            // we are in the left part...
            idVec[1-id][ibb_left_f+nbb_left].first = idVec[1-id][ibb_].first;
            idVec[1-id][ibb_left_f+nbb_left].second = ibb;
            ++nbb_left;
          }
          else {
            int_buf[2*nbb_right  ] = idVec[1-id][ibb_].first;
            int_buf[2*nbb_right+1] = ibb;
            ++nbb_right;
          }
        }
      
        for (int ibb = 0; ibb < nbb_right; ++ibb) {
          idVec[1-id][ibb_right_f+ibb].first = int_buf[2*ibb  ];
          idVec[1-id][ibb_right_f+ibb].second = int_buf[2*ibb+1];
        }
      
        // now push these new leaves onto the stack...
      
        // left...

        FOR_I6 assert(leafdata[2*ileaf+1][i] == -2);

        leafdata[2*ileaf+1][0] = child_bbminmax[2*id][0]; // jmin
        leafdata[2*ileaf+1][1] = child_bbminmax[2*id][1]; // kmin
        leafdata[2*ileaf+1][2] = -child_bbminmax[2*id][2]; assert(leafdata[2*ileaf+1][2] >= leafdata[2*ileaf+1][0]); // jmax
        leafdata[2*ileaf+1][3] = -child_bbminmax[2*id][3]; assert(leafdata[2*ileaf+1][3] >= leafdata[2*ileaf+1][1]); // kmax
        
        if (nbb_left == 1) {
          // this is the terminal leaf...
          leafdata[2*ileaf+1][4] = leafdata[2*ileaf+1][5] = mybbVec[5*ibb_left_f]; // index
        }
        else {
          // put the leaf on the stack...
          assert(nbb_left > 1);
          leafdata[2*ileaf+1][4] = ibb_left_f;
          leafdata[2*ileaf+1][5] = ibb_left_l;
          stack[istack++] = 2*ileaf+1;
        }
      
        // right...
      
        FOR_I6 assert(leafdata[2*ileaf+2][i] == -2);
      
        leafdata[2*ileaf+2][0] = child_bbminmax[2*id+1][0]; // jmin
        leafdata[2*ileaf+2][1] = child_bbminmax[2*id+1][1]; // kmin
        leafdata[2*ileaf+2][2] = -child_bbminmax[2*id+1][2]; assert(leafdata[2*ileaf+2][2] >= leafdata[2*ileaf+2][0]); // jmax
        leafdata[2*ileaf+2][3] = -child_bbminmax[2*id+1][3]; assert(leafdata[2*ileaf+2][3] >= leafdata[2*ileaf+2][1]); // kmax
      
        if (nbb_right == 1) {
          // this is the terminal leaf...
          leafdata[2*ileaf+2][4] = leafdata[2*ileaf+2][5] = mybbVec[5*ibb_right_l]; // index
        }
        else {
          // put the leaf on the stack...
          assert(nbb_right > 1);
          leafdata[2*ileaf+2][4] = ibb_right_f;
          leafdata[2*ileaf+2][5] = ibb_right_l;
          stack[istack++] = 2*ileaf+2;
        }
      
      }

      delete[] flag;
      delete[] int_buf;
    
      if (mpi_rank == 0) cout << " > done serial: " << MPI_Wtime()-wtime << endl;
    
    }
    
  }
  
};

#endif
