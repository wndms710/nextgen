#ifndef _SUB_SURFACE_HPP_
#define _SUB_SURFACE_HPP_

#include "Adt.hpp"
#include "Prcomm.hpp"

class SubSurface {
public:

  // TODO: uint8 would be better here, but infrastructure does not support it. E.g. 
  // used in CtiCanvas...
  int8 *ist_global_and_bits; // st_global-and-bits-of-st -- back reference to full surface...
  int8 *isp_global_and_bits; // sp_global-and-bits-of-sp -- back reference to full surface...
  
  int8 nsp_global,nst_global;

  double bbmin[3];
  double bbmax[3];

  int nsp,nst,nse;
  double (*xp)[3]; // local node copy - note that periodicity is no longer important
  int (*spost)[3]; // local sp-of-st
  int (*seost)[3]; // local se-of-st
  int *znost;

  Adt<double> * adt;

  int *sp_flag;
  int *st_flag;
  int *se_flag;

  double *sp_buf; 

  bool b_prcomms;
  vector<Prcomm> spPrcommVec; 
  vector<Prcomm> stPrcommVec; 
  map<const void*,MpiRequestStuff*> mpiRequestMap;

  SubSurface() {

    nsp_global = nst_global = 0;
    nsp = nst = nse = 0;

    xp = NULL;
    ist_global_and_bits = NULL;
    isp_global_and_bits = NULL;
    spost = NULL;
    seost = NULL;
    znost = NULL;
    
    adt = NULL;

    sp_flag = NULL; // this is used for fast node matching
    st_flag = NULL; 
    se_flag = NULL;

    sp_buf = NULL;

    b_prcomms = false;

  }

  ~SubSurface() {

    clear();
    
  }
  
  void setIstGlobalAndBits(const int ist,const int bits,const int ist_ss) {
    assert(ist_global_and_bits);
    assert((ist >= 0)&&(ist < TWO_BILLION));
    assert((bits >= 0)&&(bits < (1<<6)));
    assert((ist_ss >= 0)&&(ist_ss < nst));
    ist_global_and_bits[ist_ss] = (int8(bits)<<(52)) | ist;  // "room for growth" - CI Oct 2017
    // check...
    int ist_check,bits_check;
    getIstGlobalAndBits(ist_check,bits_check,ist_ss);
    assert(ist_check == ist);
    assert(bits_check == bits);
  }

  void getIstGlobalAndBits(int& ist,int& bits,const int ist_ss) const {
    assert(ist_global_and_bits);
    assert((ist_ss >= 0)&&(ist_ss < nst));
    bits = (ist_global_and_bits[ist_ss]>>(52));
    assert((bits >= 0)&&(bits < (1<<6)));
    ist = (ist_global_and_bits[ist_ss]&MASK_52BITS);
  }

  void setIspGlobalAndBits(const int isp,const int bits,const int isp_ss) {
    assert(isp_global_and_bits);
    assert((isp >= 0)&&(isp < TWO_BILLION));
    assert((bits >= 0)&&(bits < (1<<6)));
    assert((isp_ss >= 0)&&(isp_ss < nsp));
    isp_global_and_bits[isp_ss] = (int8(bits)<<(52)) | isp;  // "room for growth" - CI Oct 2017
    // check...
    int isp_check,bits_check;
    getIspGlobalAndBits(isp_check,bits_check,isp_ss);
    assert(isp_check == isp);
    assert(bits_check == bits);
  }

  void getIspGlobalAndBits(int& isp,int& bits,const int isp_ss) const {
    assert(isp_global_and_bits);
    assert((isp_ss >= 0)&&(isp_ss < nsp));
    bits = (isp_global_and_bits[isp_ss]>>(52));
    assert((bits >= 0)&&(bits < (1<<6)));
    isp = (isp_global_and_bits[isp_ss]&MASK_52BITS);
    assert((isp >= 0)&&(isp < TWO_BILLION));
  }
  
  void clear() {

    FOR_I3 bbmin[i] = HUGE_VAL;
    FOR_I3 bbmax[i] = -HUGE_VAL;
    
    nsp_global = nst_global = 0;
    nsp = nst = nse = 0;

    DELETE(xp);
    DELETE(ist_global_and_bits);
    DELETE(isp_global_and_bits);
    DELETE(spost);
    DELETE(seost);
    DELETE(znost);
    
    if (adt != NULL) {
      delete adt;
      adt = NULL;
    }

    DELETE(sp_flag);
    DELETE(st_flag);
    DELETE(se_flag);

    DELETE(sp_buf);

    b_prcomms = false;
    spPrcommVec.clear();
    stPrcommVec.clear();

  }

  void writeTecplot() {

    // not parallel...
    
    char filename[128];
    sprintf(filename,"subSurface.%06d.dat",mpi_rank);
    
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename);
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    
    fprintf(fp,"ZONE T=\"full surface\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);
    
    // data...
    for (int isp = 0; isp < nsp; ++isp)
      fprintf(fp,"%lf %lf %lf\n",xp[isp][0],xp[isp][1],xp[isp][2]);
    
    for (int ist = 0; ist < nst; ++ist)
      fprintf(fp,"%d %d %d\n",
	      spost[ist][0]+1,
	      spost[ist][1]+1,
	      spost[ist][2]+1);
    
    fclose(fp);
    
  }

  // updateSpData(int * s and int (*s)[3]...
  // updateStData(int * s and int (*s)[3]...

#define T int
#define PACK_BUF pack_buf_int
#define UNPACK_BUF unpack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 32124
#include "updateSubSurfaceDN.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateSpData(double * s and double (*s)[3]...
  // updateStData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define UNPACK_BUF unpack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 32125
#include "updateSubSurfaceDN.hpp"
#include "updateSubSurfaceDN3.hpp"
#undef T
#undef PACK_BUF
#undef UNPACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  void buildPrcomms() {

    assert(b_prcomms == false);
    b_prcomms = true; // use the bool because the *PrcommVecs could be empty after build...

    if (mpi_rank == 0)
      cout << "SubSurface::buildPrcomms()" << endl;

    // spPrcommVec...
    {
      assert(spPrcommVec.empty());

      int * spora_s = NULL; // s == striped
      //MiscUtils::buildUniformXora(spora_s,nsp_global);
      MiscUtils::calcThresholdDist(spora_s,nsp_global,mpi_size,DIST_THRESHOLD);
      const int nsp_s = spora_s[mpi_rank+1] - spora_s[mpi_rank];

      int* send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      int* send_disp = new int[mpi_size];
      int* send_buf_int = NULL;
      for (int iter = 0; iter < 2; ++iter) {

        for (int isp = 0; isp < nsp; ++isp) {
          int isp_global,bits; getIspGlobalAndBits(isp_global,bits,isp);
          const int rank = MiscUtils::getRankInXora(isp_global,spora_s);
          if (iter == 0) {
            send_count[rank] += 3; // rank/bits,isp,isp
          }
          else {
            send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(mpi_rank,bits);
            send_buf_int[send_disp[rank]++] = isp_global;
            send_buf_int[send_disp[rank]++] = isp;
          }
        }

        // rewind...
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        if (iter == 0) {
          int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          send_buf_int = new int[send_count_sum];
        }

      }

      // recv-side stuff...
      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int * recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int* recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,
          mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;

      // exchange rank/bits,isp's sharing an isp...

      set<uint8>* rbiosp = new set<uint8>[nsp_s]; // rank/bits/isp
      for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
        int rank,bits; BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv]);
        assert((rank >= 0)&&(rank < mpi_size));
        assert((bits >= 0)&&(bits < (1<<6)));
        const int isp_global = recv_buf_int[irecv+1]; 
        assert((isp_global >= spora_s[mpi_rank])&&(isp_global < spora_s[mpi_rank+1]));
        const int isp_s = isp_global-spora_s[mpi_rank];
        const int isp = recv_buf_int[irecv+2];
        const uint8 rbi = BitUtils::packRankBitsIndex(rank,bits,isp);
        assert(rbiosp[isp_s].find(rbi) == rbiosp[isp_s].end());
        rbiosp[isp_s].insert(rbi);
      }
      delete[] spora_s;
      delete[] recv_buf_int; recv_buf_int = NULL;

      FOR_RANK send_count[rank] = 0;
      for (int iter = 0; iter < 2; ++iter) {
        for (int isp_s = 0; isp_s < nsp_s; ++isp_s) {
          const int cnt = rbiosp[isp_s].size();
          assert(cnt >= 0); 
          for (set<uint8>::iterator it = rbiosp[isp_s].begin(); it != rbiosp[isp_s].end(); ++it) {
            int rank,bits,isp; BitUtils::unpackRankBitsIndex(rank,bits,isp,*it);
            assert((rank >= 0)&&(rank < mpi_size));
            if (iter == 0) {
              send_count[rank] += 2*cnt; // cnt,isp then rank_nbr/bits_nbr,isp_nbr for all nbrs sharing isp
            }
            else {
              send_buf_int[send_disp[rank]++] = cnt-1;
              send_buf_int[send_disp[rank]++] = isp;
              for (set<uint8>::iterator it2 = rbiosp[isp_s].begin(); it2 != rbiosp[isp_s].end(); ++it2) {
                if (*it2 == *it)
                  continue;
                int rank2,bits2,isp2; BitUtils::unpackRankBitsIndex(rank2,bits2,isp2,*it2);
                assert((rank2 >= 0)&&(rank2 < mpi_size));
                const int add_bits = BitUtils::addPeriodicBits(bits,BitUtils::flipPeriodicBits(bits2));
                send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(rank2,add_bits);
                send_buf_int[send_disp[rank]++] = isp2;
              }
            }
          }
        }

        // rewind...
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        if (iter == 0) {
          int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
        }
      }
      delete[] rbiosp;

      // recv-side stuff...
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,
          mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;
      delete[] recv_count;
      delete[] recv_disp;
      delete[] send_count;
      delete[] send_disp;

      vector<pair<uint8,int> > rbi_index_pair_vec;
      {
        int irecv = 0;
        while (irecv < recv_count_sum) {
          const int cnt = recv_buf_int[irecv++];
          const int isp = recv_buf_int[irecv++];
          for (int ii = 0; ii < cnt; ++ii) {
            int rank2,bits2; BitUtils::unpackRankBits(rank2,bits2,recv_buf_int[irecv++]);
            assert((rank2 >= 0)&&(rank2 < mpi_size));
            const int isp2 = recv_buf_int[irecv++];
            rbi_index_pair_vec.push_back(pair<uint8,int>(BitUtils::packRankBitsIndex(rank2,bits2,isp2),isp));
          }
        }
      }
      delete[] recv_buf_int;

      sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

      assert(spPrcommVec.empty());

      uint8 * send_buf = new uint8[rbi_index_pair_vec.size()];
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
        send_buf[ii] = rbi_index_pair_vec[ii].first;
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
        //cout << bits << " " << rank << " " << mpi_rank << " " <<  endl;
        //assert(bits||(rank != mpi_rank));
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) {
            prcomm->unpack_size = ii-prcomm->unpack_offset;
            prcomm->unpackVec.resize(prcomm->unpack_size);
            for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
              prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
            }
          }
          rank_current = rank;
          bits_current = -1;
          spPrcommVec.push_back(Prcomm());
          prcomm = &spPrcommVec.back();
          prcomm->rank = rank;
          prcomm->unpack_offset = ii;
        }
        else {
          assert(rank_current == rank);
          assert(prcomm);
        }
        if (bits > bits_current) {
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
      }
      // we just finished. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
        prcomm->unpackVec.resize(prcomm->unpack_size);
        for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
          prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
        }
      }

      // finally, we need to send/recv the indices and bits to the pack side
      // and build the packVecs. 

      MPI_Request * sendRequestArray = new MPI_Request[spPrcommVec.size()];
      MPI_Request * recvRequestArray = new MPI_Request[spPrcommVec.size()];
      for (int ii = 0,ii_end=spPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(&(spPrcommVec[ii].pack_size),1,MPI_INT,spPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

        // and the send...
        MPI_Issend(&(spPrcommVec[ii].unpack_size),1,MPI_INT,spPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(spPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(spPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      // now send from the unpack side to the pack side...

      int pack_size = 0;
      for (int ii = 0,ii_end=spPrcommVec.size(); ii < ii_end; ++ii)
        pack_size += spPrcommVec[ii].pack_size;

      uint8 * recv_rbi = new uint8[pack_size];
      pack_size = 0;
      for (int ii = 0,ii_end=spPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(recv_rbi+pack_size,spPrcommVec[ii].pack_size,MPI_UINT8,
            spPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
        pack_size += spPrcommVec[ii].pack_size;

        // and the send...
        MPI_Issend(send_buf+spPrcommVec[ii].unpack_offset,spPrcommVec[ii].unpack_size,MPI_UINT8,
            spPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(spPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(spPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      delete[] send_buf;
      delete[] sendRequestArray;
      delete[] recvRequestArray;

      // now build the packVec (and periodicity in the future)...

      double R[9], t[3];
      pack_size = 0;
      for (int ii = 0,ii_end=spPrcommVec.size(); ii < ii_end; ++ii) {
        assert(spPrcommVec[ii].packVec.empty());
        assert(spPrcommVec[ii].pack_size == spPrcommVec[ii].unpack_size);
        spPrcommVec[ii].packVec.resize(spPrcommVec[ii].pack_size);
        CvPrcomm::Transform * transform = NULL;
        int bits_current = -1;
        for (int i = 0; i < spPrcommVec[ii].pack_size; ++i) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi[pack_size+i]);
          assert(rank == mpi_rank);
          if (bits > bits_current) {
            // bits are about to change, so complete any translate and rotate...
            if (transform) transform->end = i;
            // look for new rotations/translations...
            const bool has_R = PeriodicData::getPeriodicR(R,bits);
            const bool has_t = PeriodicData::getPeriodicT(t,bits);
            spPrcommVec[ii].transformVec.push_back(Prcomm::Transform(has_R,R,has_t,t,bits,i));
            transform = &(spPrcommVec[ii].transformVec.back());
            bits_current = bits;
          }
          else {
            assert(bits_current == bits);
          }
          assert((index >= 0)&&(index < nsp));
          spPrcommVec[ii].packVec[i] = index;
        }
        if (transform) transform->end = spPrcommVec[ii].pack_size;
        pack_size += spPrcommVec[ii].pack_size;
      }
      delete[] recv_rbi;

      // =====================================================
      // test that the node coordinates are exactly matched...
      // =====================================================

      updateSpDataCheck(xp);
    }

    // stPrcommVec...
    {
      assert(stPrcommVec.empty());

      int * stora_s = NULL; // s == striped
      //MiscUtils::buildUniformXora(stora_s,nst_global);
      MiscUtils::calcThresholdDist(stora_s,nst_global,mpi_size,DIST_THRESHOLD);
      const int nst_s = stora_s[mpi_rank+1] - stora_s[mpi_rank];

      int* send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      int* send_disp = new int[mpi_size];
      int* send_buf_int = NULL;
      for (int iter = 0; iter < 2; ++iter) {

        for (int ist = 0; ist < nst; ++ist) {
          int ist_global,bits; getIstGlobalAndBits(ist_global,bits,ist);
          const int rank = MiscUtils::getRankInXora(ist_global,stora_s);
          if (iter == 0) {
            send_count[rank] += 3; // rank/bits,ist,ist
          }
          else {
            send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(mpi_rank,bits);
            send_buf_int[send_disp[rank]++] = ist_global;
            send_buf_int[send_disp[rank]++] = ist;
          }
        }

        // rewind...
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        if (iter == 0) {
          int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          send_buf_int = new int[send_count_sum];
        }

      }

      // recv-side stuff...
      int * recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int * recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int* recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,
          mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;

      // exchange rank/bits,ist's sharing an ist...

      set<uint8>* rbiost = new set<uint8>[nst_s]; // rank/bits/ist
      for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
        int rank,bits; BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv]);
        assert((rank >= 0)&&(rank < mpi_size));
        assert((bits >= 0)&&(bits < (1<<6)));
        const int ist_global = recv_buf_int[irecv+1]; 
        assert((ist_global >= stora_s[mpi_rank])&&(ist_global < stora_s[mpi_rank+1]));
        const int ist_s = ist_global-stora_s[mpi_rank];
        const int ist = recv_buf_int[irecv+2];
        const uint8 rbi = BitUtils::packRankBitsIndex(rank,bits,ist);
        assert(rbiost[ist_s].find(rbi) == rbiost[ist_s].end());
        rbiost[ist_s].insert(rbi);
      }
      delete[] stora_s;
      delete[] recv_buf_int; recv_buf_int = NULL;

      FOR_RANK send_count[rank] = 0;
      for (int iter = 0; iter < 2; ++iter) {
        for (int ist_s = 0; ist_s < nst_s; ++ist_s) {
          const int cnt = rbiost[ist_s].size();
          assert(cnt >= 0); 
          for (set<uint8>::iterator it = rbiost[ist_s].begin(); it != rbiost[ist_s].end(); ++it) {
            int rank,bits,ist; BitUtils::unpackRankBitsIndex(rank,bits,ist,*it);
            assert((rank >= 0)&&(rank < mpi_size));
            if (iter == 0) {
              send_count[rank] += 2*cnt; // cnt,ist then rank_nbr/bits_nbr,ist_nbr for all nbrs sharing ist
            }
            else {
              send_buf_int[send_disp[rank]++] = cnt-1;
              send_buf_int[send_disp[rank]++] = ist;
              for (set<uint8>::iterator it2 = rbiost[ist_s].begin(); it2 != rbiost[ist_s].end(); ++it2) {
                if (*it2 == *it)
                  continue;
                int rank2,bits2,ist2; BitUtils::unpackRankBitsIndex(rank2,bits2,ist2,*it2);
                assert((rank2 >= 0)&&(rank2 < mpi_size));
                const int add_bits = BitUtils::addPeriodicBits(bits,BitUtils::flipPeriodicBits(bits2));
                send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(rank2,add_bits);
                send_buf_int[send_disp[rank]++] = ist2;
              }
            }
          }
        }

        // rewind...
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        if (iter == 0) {
          int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
        }
      }
      delete[] rbiost;

      // recv-side stuff...
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,
          mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;
      delete[] recv_count;
      delete[] recv_disp;
      delete[] send_count;
      delete[] send_disp;

      vector<pair<uint8,int> > rbi_index_pair_vec;
      {
        int irecv = 0;
        while (irecv < recv_count_sum) {
          const int cnt = recv_buf_int[irecv++];
          const int ist = recv_buf_int[irecv++];
          for (int ii = 0; ii < cnt; ++ii) {
            int rank2,bits2; BitUtils::unpackRankBits(rank2,bits2,recv_buf_int[irecv++]);
            assert((rank2 >= 0)&&(rank2 < mpi_size));
            const int ist2 = recv_buf_int[irecv++];
            rbi_index_pair_vec.push_back(pair<uint8,int>(BitUtils::packRankBitsIndex(rank2,bits2,ist2),ist));
          }
        }
      }
      delete[] recv_buf_int;

      sort(rbi_index_pair_vec.begin(),rbi_index_pair_vec.end());

      assert(stPrcommVec.empty());

      uint8 * send_buf = new uint8[rbi_index_pair_vec.size()];
      int rank_current = -1;
      int bits_current = -1;
      Prcomm * prcomm = NULL;
      for (int ii = 0, ii_end=rbi_index_pair_vec.size(); ii < ii_end; ++ii) {
        send_buf[ii] = rbi_index_pair_vec[ii].first;
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_index_pair_vec[ii].first);
        //cout << bits << " " << rank << " " << mpi_rank << " " <<  endl;
        //assert(bits||(rank != mpi_rank));
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm) {
            prcomm->unpack_size = ii-prcomm->unpack_offset;
            prcomm->unpackVec.resize(prcomm->unpack_size);
            for (int ii2 = prcomm->unpack_offset; ii2 < ii; ++ii2) {
              prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
            }
          }
          rank_current = rank;
          bits_current = -1;
          stPrcommVec.push_back(Prcomm());
          prcomm = &stPrcommVec.back();
          prcomm->rank = rank;
          prcomm->unpack_offset = ii;
        }
        else {
          assert(rank_current == rank);
          assert(prcomm);
        }
        if (bits > bits_current) {
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
      }
      // we just finished. If there was a prcomm, then
      // complete the size...
      if (prcomm) {
        prcomm->unpack_size = rbi_index_pair_vec.size()-prcomm->unpack_offset;
        prcomm->unpackVec.resize(prcomm->unpack_size);
        for (int ii2 = prcomm->unpack_offset, ii2_end=rbi_index_pair_vec.size(); ii2 < ii2_end; ++ii2) {
          prcomm->unpackVec[ii2-prcomm->unpack_offset] = rbi_index_pair_vec[ii2].second;
        }
      }

      // finally, we need to send/recv the indices and bits to the pack side
      // and build the packVecs. 

      MPI_Request * sendRequestArray = new MPI_Request[stPrcommVec.size()];
      MPI_Request * recvRequestArray = new MPI_Request[stPrcommVec.size()];
      for (int ii = 0,ii_end=stPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(&(stPrcommVec[ii].pack_size),1,MPI_INT,stPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

        // and the send...
        MPI_Issend(&(stPrcommVec[ii].unpack_size),1,MPI_INT,stPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(stPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(stPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      // now send from the unpack side to the pack side...

      int pack_size = 0;
      for (int ii = 0,ii_end=stPrcommVec.size(); ii < ii_end; ++ii)
        pack_size += stPrcommVec[ii].pack_size;

      uint8 * recv_rbi = new uint8[pack_size];
      pack_size = 0;
      for (int ii = 0,ii_end=stPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(recv_rbi+pack_size,stPrcommVec[ii].pack_size,MPI_UINT8,
            stPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
        pack_size += stPrcommVec[ii].pack_size;

        // and the send...
        MPI_Issend(send_buf+stPrcommVec[ii].unpack_offset,stPrcommVec[ii].unpack_size,MPI_UINT8,
            stPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(stPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(stPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      delete[] send_buf;
      delete[] sendRequestArray;
      delete[] recvRequestArray;

      // now build the packVec (and periodicity in the future)...

      double R[9], t[3];
      pack_size = 0;
      for (int ii = 0,ii_end=stPrcommVec.size(); ii < ii_end; ++ii) {
        assert(stPrcommVec[ii].packVec.empty());
        assert(stPrcommVec[ii].pack_size == stPrcommVec[ii].unpack_size);
        stPrcommVec[ii].packVec.resize(stPrcommVec[ii].pack_size);
        CvPrcomm::Transform * transform = NULL;
        int bits_current = -1;
        for (int i = 0; i < stPrcommVec[ii].pack_size; ++i) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi[pack_size+i]);
          assert(rank == mpi_rank);
          if (bits > bits_current) {
            // bits are about to change, so complete any translate and rotate...
            if (transform) transform->end = i;
            // look for new rotations/translations...
            const bool has_R = PeriodicData::getPeriodicR(R,bits);
            const bool has_t = PeriodicData::getPeriodicT(t,bits);
            stPrcommVec[ii].transformVec.push_back(Prcomm::Transform(has_R,R,has_t,t,bits,i));
            transform = &(stPrcommVec[ii].transformVec.back());
            bits_current = bits;
          }
          else {
            assert(bits_current == bits);
          }
          assert((index >= 0)&&(index < nst));
          stPrcommVec[ii].packVec[i] = index;
        }
        if (transform) transform->end = stPrcommVec[ii].pack_size;
        pack_size += stPrcommVec[ii].pack_size;
      }
      delete[] recv_rbi;

      // =====================================================
      // test that the node coordinates are exactly matched...
      // =====================================================

      double (*xst)[3] = new double[nst][3];
      for (int ist = 0; ist < nst; ++ist) {
        FOR_J3 xst[ist][j] = 0.0;
        FOR_I3 {
          const int isp = spost[ist][i];
          FOR_J3 xst[ist][j] += xp[isp][j];
        }
        FOR_J3 xst[ist][j] /= 3.0;
      }
      updateStDataCheck(xst);
      delete[] xst;
    }

  }


};

#endif

