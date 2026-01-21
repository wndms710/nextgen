#ifndef UTILS_HPP
#define UTILS_HPP

#include "MiscUtils.hpp"
#include "CTI.hpp"

using namespace CTI;

inline bool icase_comp(unsigned char _a, unsigned char _b) {
  return std::tolower(_a) == std::tolower(_b);
}

inline bool equals_ci(const string& s1, const string &s2) {
  if ( s1.length() == s2.length())
    return std::equal(s1.begin(),s1.end(),s2.begin(),icase_comp);
  else
    return false;
}

// GROUP_BITS defines how many bits are used for group indexing, so max number of groups is (1<<GROUP_BITS)
// orinally this was set to 2, yielding 4 groups, but many cases arise where ngr > 4, so now using 3.
// If this becomes too large, then another strategy is required to store this info because the other
// parts of the bits will not have enough precision to store large indices.

#define GROUP_BITS 3
#define NGR_MAX (1<<(GROUP_BITS))
#define RANK_MASK 4194303

namespace BitUtils {

  // NEW stuff...
  // NEW stuff...
  // NEW stuff...

  inline void runBitTest() {
    int index = 0;
    cout << "BitUtils::runBitTest for GROUP_BITS=" << GROUP_BITS << endl;
    while (1) {
      index += rand()%1000000;
      cout << "index: " << index;
      for (int igr = 0; igr < (1<<GROUP_BITS); ++igr) {
	//assert(index < (1<<(31-GROUP_BITS)));
	int group_index = (igr<<(31-GROUP_BITS))|index;
	int igr_check = (group_index>>(31-GROUP_BITS));
	int index_check = (group_index&(~(((1<<GROUP_BITS)-1)<<(31-GROUP_BITS))));
	if ((igr_check != igr)||(index_check != index)) {
	  cout << " failed at index: " << index << " igr: " << igr << endl;
	  if (!(index < (1<<(31-GROUP_BITS))))
	    cout << "this would have been caught by index check" << endl;
	  else
	    cout << "this would NOT have been caught by index check -- must fix." << endl;
	  assert(0);
	}
	//  cout << "igr: " << igr << ", index: " << index << ", group_index: " << group_index <<
	// ", recover igr: " << (group_index>>(31-GROUP_BITS)) << ", recover index: " << (group_index&(~(((1<<GROUP_BITS)-1)<<(31-GROUP_BITS)))) << endl;
      }
      cout << " OK." << endl;
    }
  }

  // some routines use an int to

  inline int packGroupRankBits(const int ig,const int irb) { // rank-bits already together
    assert(GROUP_BITS == 3);
    assert((ig >= 0)&&(ig < (1<<GROUP_BITS)));
    assert((irb >= 0)&&(irb < (1<<(31-GROUP_BITS))));
    const int igrb = ((ig<<(31-GROUP_BITS))|irb);
    assert( (igrb>>(31-GROUP_BITS)) == ig );
    assert( (igrb&(~(((1<<GROUP_BITS)-1)<<(31-GROUP_BITS)))) == irb );
    return igrb;
  }

  inline void unpackGroupRankBits(int& ig,int& irb,const int igrb) {
    assert(GROUP_BITS == 3);
    ig = (igrb>>(31-GROUP_BITS));
    irb = (igrb&(~(((1<<GROUP_BITS)-1)<<(31-GROUP_BITS))));
    assert(packGroupRankBits(ig,irb) == igrb);
  }

  inline int packRankBits(const int rank,const int bits) {
    assert(GROUP_BITS == 3);
    // bits require 6, use the rest for rank...
    assert((bits >= 0)&&(bits < (1<<6)));
    assert((rank >= 0)&&(rank <= RANK_MASK));
    return (bits<<22) | rank;
  }

  inline void unpackRankBits(int& rank,int& bits,const int rb) {
    assert(GROUP_BITS == 3);
    bits = rb>>22; assert((bits >= 0)&&(bits < (1<<6)));
    rank = rb & RANK_MASK;
    assert(rb == packRankBits(rank,bits));
  }

  // use these for group and index...

  inline int packGroupIndex(const int ig,const int index) {
    return packGroupRankBits(ig,index);
  }

  inline void unpackGroupIndex(int& ig,int& index,const int igi) {
    unpackGroupRankBits(ig,index,igi);
  }

  // end of NEW stuff
  // end of NEW stuff
  // end of NEW stuff

  inline int8 packRankIndex(const int rank,const int index) {
    return( (int8(rank)<<32) | int8(index) );
  }

  inline void unpackRankIndex(int &rank,int &index,const int8 hash) {
    index = int( hash & ((int8(1)<<32)-1) );
    rank  = int( hash>>32 );
    assert( packRankIndex(rank,index) == hash );
  }

  inline void unpackRankBitsIndex(int& rank,int& bits,int& index,const uint8 hsh) {
    // rank is most significant, then bits, then index in first 30 bits.
    // Note this order is different than that currently used in stitch.
    // This order results in properly sorted unpack
    // ghost arrays in the maps...
    index = int( hsh & ((uint8(1)<<30)-1) );
    assert((index >= 0)&&(index < (1<<30)));
    bits = int( (hsh>>30) & ((uint8(1)<<12)-1) );
    assert((bits >= 0)&&(bits < (1<<12))); // 6 bits (for 2 pairs x 3 possible directions of periodicity) x 2
    rank = int( hsh>>(30+12) );
    assert((rank >= 0)&&(rank < (1<<22))); // allows 4M ranks
  }

  inline uint8 packRankBitsIndex(const int rank,const int bits,const int index) {
    // rank/bits/index hash...
    assert((rank >= 0)&&(rank < (1<<22))); // allows 4M ranks
    assert((bits >= 0)&&(bits < (1<<12))); // 6 bits (for 2 pairs x 3 possible directions of periodicity) x 2
    assert((index >= 0)&&(index < (1<<30)));
    uint8 hsh = int8(index);
    hsh |= (int8(bits)<<30);
    hsh |= (int8(rank)<<(30+12));
    // make sure last bit is free...
    assert(!(hsh & (uint8(1)<<63)));
    // and check that this unpacks correctly...
    int rank_check,bits_check,index_check;
    unpackRankBitsIndex(rank_check,bits_check,index_check,hsh);
    assert(rank_check == rank);
    assert(bits_check == bits);
    assert(index_check == index);
    return(hsh);
  }

  inline uint8 packRankBitsIndex(const int rank_bits,const int index) {
    int rank,bits;
    unpackRankBits(rank,bits,rank_bits);
    return packRankBitsIndex(rank,bits,index);
  }

  /*
  // part-based, surface hashing
  inline void unpackPbiHash(int& bits,int& index,const uint32_t hsh) {
    // rank is first 22 bits...
    index = ( hsh & ((1<<(32-6))-1) );
    bits  = ( hsh >> (32-6) );
  }

  inline uint32_t makePbiHash(const int bits,const int index) {
    // rank/bits/index hash...
    assert((bits >= 0)&&(bits < (1<<6)-1)); // allows 4M ranks
    assert((index >= 0)&&(index < (1<<(32-6))));
    uint32_t hsh = ( index | (bits<<(32-6)) );
    // and check that this unpacks correctly...
    int bits_check,index_check;
    unpackPbiHash(bits_check,index_check,hsh);
    assert(bits_check == bits);
    assert(index_check == index);
    return(hsh);
  }
  */

  // part-based, surface hashing: periodic bits and index...
  inline void unpackPbiHash(int& bits,int& index,const uint8 hsh) {
    // periodic bits are stored in the last 12...all others are the index
    index = ( hsh & ((uint8(1)<<(64-12))-1) );
    bits  = ( hsh >> (64-12) );
  }

  inline uint8 packPbiHash(const int bits,const int index) {
    // bits/index hash..
    assert((bits >= 0)&&(bits < (1<<12)-1));
    //assert((index >= 0)&&(index < (uint8(1)<<(64-12))));
    assert((index >= 0)&&(index < TWO_BILLION));
    uint8 hsh = ( index | (uint8(bits)<<(64-12)) );
    // and check that this unpacks correctly...
    int bits_check;
    int index_check;
    unpackPbiHash(bits_check,index_check,hsh);
    assert(bits_check == bits);
    assert(index_check == index);
    return(hsh);
  }

  // same pbi stuff, but for int8 indices...
  inline void unpackPbiHash(int& bits,int8& index,const uint8 hsh) {
    // periodic bits are stored in the last 12...all others are the index
    index = ( hsh & ((uint8(1)<<(64-12))-1) );
    bits  = ( hsh >> (64-12) );
  }

  inline uint8 packPbiHash(const int bits,const int8 index) {
    // bits/index hash..
    assert((bits >= 0)&&(bits < (1<<12)-1));
    //assert((index >= 0)&&(index < (uint8(1)<<(64-12))));
    assert((index >= 0)&&(index < int8(uint8(1)<<(64-12))));
    uint8 hsh = ( index | (uint8(bits)<<(64-12)) );
    // and check that this unpacks correctly...
    int bits_check;
    int8 index_check;
    unpackPbiHash(bits_check,index_check,hsh);
    assert(bits_check == bits);
    assert(index_check == index);
    return(hsh);
  }

  inline void unpackPartIndex(int& ip,int& index,const uint8 piv) {
    // index is first 32 bits...
    index = piv & ((uint8(1)<<32)-1);
    ip = (piv>>32);
  }

  inline uint8 makePartIndex(const int ip,const int index) {
    assert(ip >= 0);
    assert(index >= 0);
    uint8 piv = uint8(index);
    piv |= (uint8(ip)<<32);
    // and check that this unpacks correctly...
    int ip_check,index_check;
    unpackPartIndex(ip_check,index_check,piv);
    assert(ip_check == ip);
    assert(index_check == index);
    return(piv);
  }

  inline void dumpBits(const int bits,const string& message) {
    cout << message << ": ";
    for (int i = 0; i < 32; ++i)
      if (bits & (1<<i))
	cout << "1";
      else
	cout << "0";
    cout << endl;
  }

  inline void dumpBits(const int bits,const int n) {
    for (int i = 0; i < n; ++i)
      if (bits & (1<<i))
	cout << "1";
      else
	cout << "0";
    cout << endl;
  }

  // periodic bits occur in bit pairs...
  inline int flipPeriodicBits(const int bits) {
    int new_bits = 0;
    int bit_pair = 0;
    while (bits >= (1<<bit_pair)) {
      bool b0 = (bits & (1<<(bit_pair)));
      bool b1 = (bits & (1<<(bit_pair+1)));
      if (b0) {
	assert(!b1);
	new_bits |= (1<<(bit_pair+1));
      }
      else if (b1) {
	assert(!b0);
	new_bits |= (1<<(bit_pair));
      }
      // increment the bit_pair...
      bit_pair += 2;
    }
    return(new_bits);
  }
  inline int flipFirstSixPeriodicBits(const int bits) {
    int new_bits = 0;
    int bit_pair = 0;
    while ((bits >= (1<<bit_pair))&&(bit_pair < 6)) {
      bool b0 = (bits & (1<<(bit_pair)));
      bool b1 = (bits & (1<<(bit_pair+1)));
      if (b0) {
	assert(!b1);
	new_bits |= (1<<(bit_pair+1));
      }
      else if (b1) {
	assert(!b0);
	new_bits |= (1<<(bit_pair));
      }
      // increment the bit_pair...
      bit_pair += 2;
    }
    // keep latter bits
    new_bits |= (~MASK_6BITS)&bits;

    /*
    cout << "orig bits: ";
    for (int i = 11; i >= 0; --i) {
      if (bits&(1<<i)) cout << "1";
      else cout << "0";
    }
    cout << " flip bits: ";
    for (int i = 11; i >= 0; --i) {
      if (new_bits&(1<<i)) cout << "1";
      else cout << "0";
    }
    cout << endl;
    */

    return(new_bits);
  }

  inline int addPeriodicBits(const int bits0,const int bits1) {
    int new_bits = 0;
    int bit_pair = 0;
    assert(max(bits0,bits1) < (1<<6)); // make sure no double bits are set
    while (max(bits0,bits1) >= (1<<bit_pair)) {
      const bool b00 = (bits0 & (1<<(bit_pair)));
      const bool b01 = (bits0 & (1<<(bit_pair+1)));
      const bool b10 = (bits1 & (1<<(bit_pair)));
      const bool b11 = (bits1 & (1<<(bit_pair+1)));
      // there should not be both bits set in either -- this indicates a
      // problem, because these should have cancelled...
      assert(!(b00 && b01));
      assert(!(b10 && b11));
      if (b00) {
	if (b10) {
	  assert(bit_pair < 6);
	  new_bits |= (1<<(bit_pair+6)); // 3 pairs of periodicity are allowed, so shift by 6
	}
	else if (!b11) {
	  new_bits |= (1<<(bit_pair));
	}
      }
      else if (b01) {
	if (b11) {
	  assert(bit_pair < 6);
	  new_bits |= (1<<(bit_pair+7)); // 3 pairs of periodicity are allowed, so shift by 6+1
	}
	else if (!b10) {
	  new_bits |= (1<<(bit_pair+1));
	}
      }
      else if (b10) {
	new_bits |= (1<<(bit_pair));
      }
      else if (b11) {
	new_bits |= (1<<(bit_pair+1));
      }
      // increment the bit_pair...
      bit_pair += 2;
    }

    /*
      cout << "adding bits: ";
      for (int i = 11; i >= 0; --i) {
      if (bits0&(1<<i)) cout << "1";
      else cout << "0";
      }
      cout << " and ";
      for (int i = 11; i >= 0; --i) {
      if (bits1&(1<<i)) cout << "1";
      else cout << "0";
      }
      cout << " and got ";
      for (int i = 11; i >= 0; --i) {
      if (new_bits&(1<<i)) cout << "1";
      else cout << "0";
      }
      cout << endl;
    */

    //assert(new_bits < (1<<6));
    return(new_bits);
  }

  inline bool hasAnyOddBit(const int bits) {
    int bit = 1;
    while (bits >= (1<<bit)) {
      if (bits & (1<<bit))
	return true;
      bit += 2;
    }
    return false;
  }
};

// needs to persist for NGR_MAX
//#undef GROUP_BITS

inline void buildXoraNew(int8 * &Xora,const int nX_local,MPI_Comm& mpi_comm) {

  assert(Xora == NULL);

  int mpi_size;
  MPI_Comm_size(mpi_comm, &mpi_size);
  //cout << "got size: " << mpi_size << endl;

  Xora = new int8[mpi_size+1];
  int8 nX_i8 = (int8)nX_local;

  MPI_Allgather(&nX_i8,1,MPI_INT8,Xora+1,1,MPI_INT8,mpi_comm);
  Xora[0] = 0;
  for (int i = 0; i < mpi_size; ++i)
    Xora[i+1] += Xora[i];
}

inline void buildNroraInternode(int * &nrora_internode) {

  assert(nrora_internode == NULL);
  nrora_internode = new int[mpi_size_internode];
  if (mpi_rank_shared == 0) {
    MPI_Allgather(&mpi_size_shared,1,MPI_INT,nrora_internode,1,MPI_INT,mpi_comm_internode);
  }
  MPI_Bcast(nrora_internode,mpi_size_internode,MPI_INT,0,mpi_comm_shared);

}

inline void convertXoraInternodeToXora(int8 * &Xora,int * nrora_internode = NULL) {

  // this routine involves less communication if we passed in nrora_internode,
  // i.e. the rank-count of each internode_rank. This is neccessary because we
  // cannot assume that every rank has the same number of ranks/node as we do...

  int * nrora_internode_ = NULL;
  if (nrora_internode == NULL) {
    buildNroraInternode(nrora_internode_);
  }
  else {
    nrora_internode_ = nrora_internode;
  }

  if (mpi_rank_shared == 0) {
    assert(Xora != NULL);
    assert(Xora[0] == 0);
    MiscUtils::resize(Xora,mpi_size_internode+1,mpi_size+1);
  }
  else {
    assert(Xora == NULL);
    Xora = new int8[mpi_size+1];
  }
  MPI_Bcast(Xora+1,mpi_size_internode,MPI_INT8,0,mpi_comm_shared);

  int rank = mpi_size;
  for (int rank_internode = mpi_size_internode; rank_internode > 0; --rank_internode) {
    for (int ii = 0; ii < nrora_internode_[rank_internode-1]; ++ii) {
      Xora[rank--] = Xora[rank_internode];
    }
  }
  assert(rank == 0);
  Xora[0] = 0;

  if (nrora_internode == NULL)
    delete[] nrora_internode_;

}

inline void MPI_Pause(const string& message,const MPI_Comm& mpi_comm) {

  int mpi_rank;
  MPI_Comm_rank(mpi_comm, &mpi_rank);

  if (mpi_rank == 0) {
    cout << "MPI_Pause() " << message << endl;
    getchar();
  }
  MPI_Barrier(mpi_comm);

}
#endif
