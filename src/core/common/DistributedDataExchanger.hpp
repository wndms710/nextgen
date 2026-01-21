#ifndef DISTRIBUTED_DATA_EXCHANGER_HPP
#define DISTRIBUTED_DATA_EXCHANGER_HPP

#include "Common.hpp"
#include "MpiStuff.hpp"
//using namespace MpiStuff; // do NOT use the namespace -- we allow mpi_comm/rank/size to be different

// ========================================================================
// A DistributedDataExchanger can be used when the data knows where it needs to
// go or come from in terms of a single striped distribution across processors.
// Once constructed, the main commands are push() and pull(): Push implies
// pushing to the distributed data, so the distributed array is being filled. This
// operation inherently requires an action because multiple data can be pushed to
// the same location of the distributed array. The pull operation implies a pull
// from distributed data. In this case, no action is required, because the pull
// is inherently single-valued.
//
// WARNING: when using ADD_DATA, MAX_DATA, etc, be sure to initialize the
// array you are populating. The push routine does NOT do this.
// ========================================================================

#define DDE_ALLTOALL 1

class PeerNbr {
public:
  int nbr_rank ;
  int n ;
  int offset ;
  PeerNbr * next ;

public:
  PeerNbr(int _rank, int _n, int _offset) {
    nbr_rank = _rank ;
    n        = _n    ;
    offset   = _offset;
    next     = NULL ;
  }

  PeerNbr(const PeerNbr& _pr) {
    nbr_rank = _pr.nbr_rank ;
    n        = _pr.n ;
    offset   = _pr.offset ;
    next     = _pr.next ;
  }

};

template<class T>
class LinkedList {
public :
  T * root ;
  int size_  ;

  LinkedList() {
    root = NULL ;
    size_ = 0 ;
  }

public :

  void clear() {

    T * nextNode = root ;
    while ( nextNode != NULL ) {
      T * currNode = nextNode ;
      nextNode     = nextNode->next ;
      delete currNode ;
    }

    root = NULL ;
    size_ = 0;
  }

  // this is an O(N) insert ...
  void push_back(const T& node ) {
    T * myNode = new T(node) ;
    if ( root == NULL ) {
      root = myNode ;
    }
    else {
      T* nextNode = root ;
      while ( true ) {
        if ( nextNode->next == NULL )
          break;

        nextNode = nextNode->next ;
      }

      nextNode->next = myNode ;
    }
    ++size_ ;
  }


  void push_front(const T& node) {
    T * myNode = new T(node) ;
    myNode->next = root ;
    root         = myNode ;
    ++size_ ;
  }

  T* front() const {
    return root ;
  }

  int size() const {
    return size_ ;
  }


  ~LinkedList() {
    this->clear();
  }
};

//typedef vector<PeerNbr> NbrVec ;
//#define FOR_SENDNBR for (NbrVec::iterator sendIt = sendNbrVec.begin(); sendIt != sendNbrVec.end(); ++sendIt)
//#define FOR_RECVNBR for (NbrVec::iterator recvIt = recvNbrVec.begin(); recvIt != recvNbrVec.end(); ++recvIt)


typedef LinkedList<PeerNbr> NbrVec ;
#define FOR_SENDNBR  for( PeerNbr* sendIt = sendNbrVec.root ; sendIt != NULL ; sendIt = sendIt->next )
#define FOR_RECVNBR  for( PeerNbr* recvIt = recvNbrVec.root ; recvIt != NULL ; recvIt = recvIt->next )

class DistributedDataExchanger {

  friend class DistributedDataExchangerReducer;

public:

  MPI_Comm mpi_comm;
  int mpi_rank,mpi_size;

  int * send_count;
  int * send_disp;
  int send_count_sum;
  int * send_index;
  int my_nda;

  int * recv_count;
  int * recv_disp;
  int recv_count_sum;
  int * recv_index;
  int nda;

  //vector<PeerNbr> sendNbrVec ;
  //vector<PeerNbr> recvNbrVec ;
  NbrVec sendNbrVec ;
  NbrVec recvNbrVec ;


public:

  DistributedDataExchanger(const MPI_Comm& comm = MpiStuff::mpi_comm) {
    MPI_Comm_dup(comm,&mpi_comm);
    MPI_Comm_rank(comm,&mpi_rank);
    MPI_Comm_size(comm,&mpi_size);

    send_count = NULL;
    send_disp  = NULL;
    send_index = NULL;
    recv_count = NULL;
    recv_disp  = NULL;
    recv_index = NULL;

    sendNbrVec.clear() ;
    recvNbrVec.clear() ;

  }

  DistributedDataExchanger(const int * my_ida_global,const int my_nda,const int daora[],const MPI_Comm& comm = MpiStuff::mpi_comm) {
    MPI_Comm_dup(comm,&mpi_comm);
    MPI_Comm_rank(comm,&mpi_rank);
    MPI_Comm_size(comm,&mpi_size);

    init(my_ida_global,my_nda,daora);
    clearSendIndexIfPossible();
  }

  DistributedDataExchanger(const int * my_ida_global,const int my_nda,int nda_striped,const MPI_Comm& comm = MpiStuff::mpi_comm) {
    MPI_Comm_dup(comm,&mpi_comm);
    MPI_Comm_rank(comm,&mpi_rank);
    MPI_Comm_size(comm,&mpi_size);

    // note: int nda_striped should be const, but it messes with NoMpi...

    // this constructor builds the daora from the passed nda_striped, then 
    // deletes it below after it is not req'd anymore...
    int * daora = new int[mpi_size+1];
    MPI_Allgather(&nda_striped,1,MPI_INT,daora+1,1,MPI_INT,mpi_comm);
    daora[0] = 0;
    for (int i = 0; i < mpi_size; i++)
      daora[i+1] += daora[i];
    
    init(my_ida_global,my_nda,daora);
    clearSendIndexIfPossible();

    delete[] daora;
  }

  DistributedDataExchanger(const int8 * my_ida_global,const int my_nda,const int8 daora[],const MPI_Comm& comm = MpiStuff::mpi_comm) {
    MPI_Comm_dup(comm,&mpi_comm);
    MPI_Comm_rank(comm,&mpi_rank);
    MPI_Comm_size(comm,&mpi_size);

    init(my_ida_global,my_nda,daora);
    clearSendIndexIfPossible();
  }

  DistributedDataExchanger(const DistributedDataExchanger& copy) {
    UNUSED(copy);
    // do not allow a copy constructor (unless you want to write one!)...
    throw(0);
  }

  ~DistributedDataExchanger() {
    DELETE(send_count);
    DELETE(send_disp);
    DELETE(send_index);
    DELETE(recv_count);
    DELETE(recv_disp);
    DELETE(recv_index);
    MPI_Comm_free(&mpi_comm);
  }

  int getDistBufSize() const { return(recv_count_sum); }

  int getDistIndex(const int idist) const { return(recv_index[idist]); }

  void pushToDistBuf(int * recv_buf,const int * my_i0,const UpdateAction action = REPLACE_DATA) {

    int * send_buf = new int[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_i0[ida];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    doExchange<int>(send_buf,recv_buf) ;
    delete[] send_buf;
  }

  void push_one(int * i0,const UpdateAction action = REPLACE_DATA) {

    int * send_buf = new int[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = 1;
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = 0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] += 1;
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    int * recv_buf = new int[recv_count_sum];
    doExchange<int>(send_buf,recv_buf);
    delete[] send_buf;

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] = recv_buf[ir];
      }
      break;
    case ADD_DATA:
      // note that we add directly into the values already in i0...
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] += recv_buf[ir];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  void push_toCsr(int * i_v, const int * i_i, const int * my_i0) {

    int * send_buf = new int[send_count_sum];

    // the distributed array cannot have more than one
    // of the same ix_global keys for this to work ...

    for (int is =0; is < send_count_sum ; ++is)
      send_buf[is] = -1;


    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      assert ( send_buf[is] == -1) ;
      send_buf[is] = my_i0[ida];
    } //ida ..

    int * recv_buf = new int[recv_count_sum];

    doExchange<int>(send_buf,recv_buf) ;

    delete[] send_buf;

    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida = recv_index[ir];
      int dod       = i_i[ida];
      while ( dod < i_i[ida+1]) {
        if ( i_v[dod] == -1) {
          i_v[dod] = recv_buf[ir];
          break;
        }
        ++dod ;
      }
      assert ( dod != i_i[ida+1]) ;

    }//ir

    delete[] recv_buf;
  }

  void push_to_csr(int * i_i,uint8 * i_v,const uint8 * my_i0) {

    // this version uses i_i to index, then rewinds

    uint8 * send_buf = new uint8[send_count_sum];

    // the distributed array cannot have more than one
    // of the same ix_global keys for this to work ...

    for (int is =0; is < send_count_sum ; ++is)
      send_buf[is] = (uint8(1)<<63);

    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      assert (send_buf[is] == (uint8(1)<<63)) ;
      send_buf[is] = my_i0[ida];
    } //ida ..

    for (int is =0; is < send_count_sum ; ++is)
      assert(send_buf[is] != (uint8(1)<<63));

    uint8 * recv_buf = new uint8[recv_count_sum];

    doExchange<uint8>(send_buf,recv_buf) ;

    delete[] send_buf;

    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida = recv_index[ir];
      i_v[i_i[ida]++] = recv_buf[ir];
    }

    delete[] recv_buf;

    // rewind i_i...

    for (int ida = nda-1; ida > 0; --ida)
      i_i[ida] = i_i[ida-1];
    i_i[0] = 0;

  }

  void push(int * i0,const int * my_i0,const UpdateAction action = REPLACE_DATA) {

    int * send_buf = new int[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_i0[ida];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = 0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] += my_i0[ida];
      }
      break;
    case MIN_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = BIG_INT;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = min(my_i0[ida],send_buf[is]);
      }
      break;
    default:
      CERR("unsupported action in push: " << action);
    }

    int * recv_buf = new int[recv_count_sum];
    doExchange<int>(send_buf,recv_buf) ;
    delete[] send_buf;

    // note that for all these cases (and particularly add and min/max) we do NOT
    // initialize the data. i.e. we add/min/max against what is already there. The
    // calling routine has to initialize this buffer appropriately. It is not done
    // in here because in some cases we want the behaviour of operating additionally
    // on the data that is already there (e.g. averaging operator).

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] = recv_buf[ir];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] += recv_buf[ir];
      }
      break;
    case MIN_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] = min(i0[ida],recv_buf[ir]);
      }
      break;
    default:
      CERR("unsupported action in push: " << action);
    }

    delete[] recv_buf;
  }

  void push(uint2 * i0,const uint2 * my_i0,const UpdateAction action = REPLACE_DATA) {
    
    uint2 * send_buf = new uint2[send_count_sum];
    
    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_i0[ida];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = uint2(0);
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] += my_i0[ida];
      }
      break;
    case MIN_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = USHRT_MAX;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = min(my_i0[ida],send_buf[is]);
      }
      break;
    default:
      CERR("unsupported action in push: " << action);
    }
    
    uint2 * recv_buf = new uint2[recv_count_sum];
    doExchange<uint2>(send_buf,recv_buf) ; 
    delete[] send_buf;

    // note that for all these cases (and particularly add and min/max) we do NOT 
    // initialize the data. i.e. we add/min/max against what is already there. The
    // calling routine has to initialize this buffer appropriately. It is not done
    // in here because in some cases we want the behaviour of operating additionally 
    // on the data that is already there (e.g. averaging operator).

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] = recv_buf[ir];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] += recv_buf[ir];
      }
      break;
    case MIN_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	i0[ida] = min(i0[ida],recv_buf[ir]);
      }
      break;
    default:
      CERR("unsupported action in push: " << action);
    }
    
    delete[] recv_buf;
  }

  void push(int8 * d0,const int8 * my_d0,const UpdateAction action = REPLACE_DATA) {

    int8 * send_buf = new int8[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_d0[ida];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = int8(0);
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] += my_d0[ida];
      }
      break;
    case MAX_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = -(int8(1)<<62); // need a large negative number here
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = max(send_buf[is],my_d0[ida]);
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    int8 * recv_buf = new int8[recv_count_sum];
    doExchange<int8>(send_buf,recv_buf);
    delete[] send_buf;

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] = recv_buf[ir];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] += recv_buf[ir];
      }
      break;
    case MAX_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] = max(d0[ida],recv_buf[ir]);
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  void push(uint8 * d0,const uint8 * my_d0,const UpdateAction action = REPLACE_DATA) {

    uint8 * send_buf = new uint8[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_d0[ida];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    uint8 * recv_buf = new uint8[recv_count_sum];
    doExchange<uint8>(send_buf,recv_buf);
    delete[] send_buf;

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] = recv_buf[ir];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  void push(double * d0,const double * my_d0,const UpdateAction action = REPLACE_DATA) {

    double * send_buf = new double[send_count_sum];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = my_d0[ida];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = 0.0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] += my_d0[ida];
      }
      break;
    case MAX_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	send_buf[is] = -1.0E+20;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	send_buf[is] = max(send_buf[is],my_d0[ida]);
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    double * recv_buf = new double[recv_count_sum];
    doExchange<double>(send_buf,recv_buf);
    delete[] send_buf;

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] = recv_buf[ir];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] += recv_buf[ir];
      }
      break;
    case MAX_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	d0[ida] = max(d0[ida],recv_buf[ir]);
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  void pushToDistBuf(double (*recv_buf)[3],const double (*my_d0)[3],const UpdateAction action = REPLACE_DATA) {

    double * send_buf = new double[send_count_sum*3];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I3 send_buf[3*is+i] = my_d0[ida][i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }

    double * recv_buf_tmp = (double*) recv_buf;
    doExchange<double>(send_buf,recv_buf_tmp) ;
    delete[] send_buf;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }
  }

  template<typename T>
  void push( T (*d0)[3],const T (*my_d0)[3],const UpdateAction action = REPLACE_DATA) {

    T * send_buf = new T[send_count_sum*3];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I3 send_buf[3*is+i] = my_d0[ida][i];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	FOR_I3 send_buf[3*is+i] = 0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I3 send_buf[3*is+i] += my_d0[ida][i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    FOR_SENDNBR {
      sendIt->n *= 3;
      sendIt->offset *= 3;
    }

    FOR_RECVNBR {
      recvIt->n *= 3;
      recvIt->offset *= 3;
    }

    T* recv_buf = new T[recv_count_sum*3];
    doExchange<T>(send_buf,recv_buf) ;

    delete[] send_buf;

    FOR_SENDNBR {
      sendIt->n /= 3;
      sendIt->offset /= 3;
    }

    FOR_RECVNBR {
      recvIt->n /= 3;
      recvIt->offset /= 3;
    }

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I3 d0[ida][i] = recv_buf[3*ir+i];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I3 d0[ida][i] += recv_buf[3*ir+i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  template<typename T>
  void push( T (*d0)[2],const T (*my_d0)[2],const UpdateAction action = REPLACE_DATA) {

    T * send_buf = new T[send_count_sum*2];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I2 send_buf[2*is+i] = my_d0[ida][i];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	FOR_I2 send_buf[2*is+i] = 0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I2 send_buf[2*is+i] += my_d0[ida][i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    FOR_SENDNBR {
      sendIt->n *= 2;
      sendIt->offset *= 2;
    }

    FOR_RECVNBR {
      recvIt->n *= 2;
      recvIt->offset *= 2;
    }

    T* recv_buf = new T[recv_count_sum*2];
    doExchange<T>(send_buf,recv_buf) ;

    delete[] send_buf;

    FOR_SENDNBR {
      sendIt->n /= 2;
      sendIt->offset /= 2;
    }

    FOR_RECVNBR {
      recvIt->n /= 2;
      recvIt->offset /= 2;
    }

    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I2 d0[ida][i] = recv_buf[2*ir+i];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I2 d0[ida][i] += recv_buf[2*ir+i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  void push(double (*d0)[3],const double (*my_d0)[3],const UpdateAction action = REPLACE_DATA) {

    double * send_buf = new double[send_count_sum*3];

    switch (action) {
    case REPLACE_DATA:
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I3 send_buf[3*is+i] = my_d0[ida][i];
      }
      break;
    case ADD_DATA:
      for (int is = 0; is < send_count_sum; ++is)
	FOR_I3 send_buf[3*is+i] = 0.0;
      for (int ida = 0; ida < my_nda; ++ida) {
	const int is = send_index[ida];
	FOR_I3 send_buf[3*is+i] += my_d0[ida][i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }


    double * recv_buf = new double[recv_count_sum*3];
    doExchange<double>(send_buf,recv_buf) ;
    delete[] send_buf;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }


    switch (action) {
    case REPLACE_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I3 d0[ida][i] = recv_buf[3*ir+i];
      }
      break;
    case ADD_DATA:
      for (int ir = 0; ir < recv_count_sum; ++ir) {
	const int ida = recv_index[ir];
	FOR_I3 d0[ida][i] += recv_buf[3*ir+i];
      }
      break;
    default:
      if (mpi_rank == 0)
	cerr << ERRSTART << "Error: unsupported action in push: " << action << ERREND << endl;
      throw(0);
    }

    delete[] recv_buf;
  }

  template<typename T>
  void push_csr(T* i_v, const T* i_i, const T* my_i_i, const T* my_i_v, const int8* ix_global,
                const int8 * daora) {

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>();

    // this is the reverse routine of pull_csr ...
    int * my_recv_count = new int[mpi_size];
    int * my_recv_disp  = new int[mpi_size];
    int * my_send_count = new int[mpi_size];
    int * my_send_disp  = new int[mpi_size];

    FOR_RANK my_send_count[rank] = 0;
    int * send_buf = new int[send_count_sum];


    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      send_buf[is]  = int(my_i_i[ida+1]-my_i_i[ida])+2; // + bookkeeping info...
    }

    FOR_SENDNBR {
      const int offset   = sendIt->offset ;
      const int count    = sendIt->n ;
      const int nbr_rank = sendIt->nbr_rank ;

      for (int i = 0 ; i < count ; ++i)
        my_send_count[nbr_rank] += send_buf[i+offset];
    }

    // build the recv side.
    MPI_Alltoall(my_send_count,1,MPI_INT,my_recv_count,1,MPI_INT,mpi_comm);

    int my_send_count_sum = 0 ;
    int my_recv_count_sum = 0 ;

    my_send_disp[0]       = 0 ;
    my_recv_disp[0]       = 0 ;
    for (int rank = 1; rank < mpi_size ; ++rank) {
      my_send_disp[rank] = my_send_disp[rank-1] + my_send_count[rank-1];
      my_send_count_sum += my_send_count[rank-1];

      my_recv_disp[rank] = my_recv_disp[rank-1] + my_recv_count[rank-1];
      my_recv_count_sum += my_recv_count[rank-1];
    }

    my_send_count_sum += my_send_count[mpi_size-1] ;
    my_recv_count_sum += my_recv_count[mpi_size-1] ;

    // let send_buf store the offset for the send..
    if ( my_send_count_sum > 0 ) {
      send_buf[my_nda-1] = my_send_count_sum - send_buf[my_nda-1];
      for (int is = my_nda-2; is >= 0; --is) {
        send_buf[is] = send_buf[is+1] - send_buf[is];
      }
      assert ( send_buf[0] == 0);
    }

    T* send_buf_T = new T[my_send_count_sum];

    // recall the send side is distributed data side...
    int count =0;
    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      int ii       = send_buf[is];
      // global index, the counts, then the ranks...
      T this_ida_global    = ix_global[ida];
      send_buf_T[ii++]     = -(this_ida_global)-1; ++count;
      send_buf_T[ii++]     = my_i_i[ida+1] - my_i_i[ida];  ++count;

      for (T dod = my_i_i[ida]; dod != my_i_i[ida+1]; ++dod) {
        send_buf_T[ii++]  = my_i_v[dod];
        ++count;
      }
    }//is

    assert ( count == my_send_count_sum ) ;
    DELETE(send_buf);

    T* recv_buf_T = new T[my_recv_count_sum];

    MPI_Alltoallv(send_buf_T,my_send_count,my_send_disp,MPI_T,
                  recv_buf_T,my_recv_count,my_recv_disp,MPI_T,mpi_comm);

    DELETE(send_buf_T);

    // now unpack on the other side...
    {
      int ii = 0 ;
      int ir = 0 ;
      while ( ii < my_recv_count_sum ) {

        const int ida       = recv_index[ir] ;
        T this_ida_global   = recv_buf_T[ii++]; assert( this_ida_global < 0 );
        this_ida_global     = -this_ida_global-1;
        assert ( ida+daora[mpi_rank] == this_ida_global);

        int nod = int(recv_buf_T[ii++]);

        // if you find a -1 entry, then we can start copying in.
        int ioi = i_i[ida];
        while ( i_v[ioi] != -1)
          ++ioi ;

        assert( i_i[ida+1]-ioi >= nod) ; // make sure theres enough space.

        for (int jj = 0 ; jj < nod ; ++jj)
          i_v[ioi+jj] = recv_buf_T[ii++];

        ++ir; // next
      }
    }

    DELETE(recv_buf_T);
    DELETE(my_recv_count);
    DELETE(my_recv_disp);
    DELETE(my_send_count);
    DELETE(my_send_disp);
  }

  void pushToCsr(int * i_v, const int * i_i, const int * my_i0) {

    int * send_buf = new int[send_count_sum];

    // the distributed array cannot have more than one
    // of the same ix_global keys for this to work ...

    for (int is =0; is < send_count_sum ; ++is)
      send_buf[is] = -1;


    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      assert ( send_buf[is] == -1) ;
      send_buf[is] = my_i0[ida];
    } //ida ..

    int * recv_buf = new int[recv_count_sum];

    doExchange<int>(send_buf,recv_buf) ;

    delete[] send_buf;

    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida = recv_index[ir];
      int dod       = i_i[ida];
      while ( dod < i_i[ida+1]) {
        if ( i_v[dod] == -1) {
          i_v[dod] = recv_buf[ir];
          break;
        }
        ++dod ;
      }
      assert ( dod != i_i[ida+1]) ;

    }//ir

    delete[] recv_buf;
  }

#ifdef WHY_ARE_BOTH_THESE_TYPE_T

  template<typename T>
  void pull_csr(T* my_i_v, const T* my_i_i, const T* i_i, const T* i_v,const int8* ix_global,
                const int8 * daora) {

    // XXX .. i need to build a map of my local indices to the global ones...
    // there may be a more efficient way of doing this ..
    /*
    std::map<int8,int> ixoxg ;
    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int8 this_ida_global = ix_global[ida];
      ixoxg[this_ida_global]     = ida ;
    }
    */

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>();

    int * my_recv_count = new int[mpi_size];
    int * my_recv_disp  = new int[mpi_size];

    int * my_send_count = new int[mpi_size];
    int * my_send_disp  = new int[mpi_size];

    FOR_RANK my_recv_count[rank] = 0;


    int* recv_buf = new int[recv_count_sum];
    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida   = recv_index[ir];
      recv_buf[ir]    = int(i_i[ida+1]-i_i[ida])+2; // take some bookkeeping info too...
    }

    FOR_RECVNBR {

      const int offset   = recvIt->offset ;
      const int count    = recvIt->n ;
      const int nbr_rank = recvIt->nbr_rank ;

      for (int i =0; i < count ; ++i) {
        my_recv_count[nbr_rank] += recv_buf[i+offset];
      }
    }//recvnbr

    // build the send side ...
    MPI_Alltoall(my_recv_count,1,MPI_INT,my_send_count,1,MPI_INT,mpi_comm);

    // build the displacements ..
    int my_recv_count_sum = 0;
    int my_send_count_sum = 0;

    my_recv_disp[0] = 0;
    my_send_disp[0] = 0;
    for (int rank =1 ; rank < mpi_size ; ++rank) {
      my_recv_disp[rank] = my_recv_disp[rank-1] + my_recv_count[rank-1];
      my_recv_count_sum += my_recv_count[rank-1];

      my_send_disp[rank] = my_send_disp[rank-1] + my_send_count[rank-1];
      my_send_count_sum += my_send_count[rank-1];
    }

    my_recv_count_sum += my_recv_count[mpi_size-1];
    my_send_count_sum += my_send_count[mpi_size-1];

    DELETE(recv_buf);

    // pack ...
    T* recv_buf_T = new T[my_recv_count_sum];

    int ii = 0;
    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida = recv_index[ir];

      // global index, the count, then the ranks ..
      T this_ida_global   = T( ida + T(daora[mpi_rank]));
      recv_buf_T[ii++] = -(this_ida_global)-1;
      recv_buf_T[ii++] = (i_i[ida+1]-i_i[ida]) ;

      for (int dod = i_i[ida]; dod != i_i[ida+1]; ++dod) {
        recv_buf_T[ii++] = i_v[dod];
      }
    }//ir

    assert ( ii == my_recv_count_sum ) ;

    T* send_buf_T = new T[my_send_count_sum];

    // do the exchange
    MPI_Alltoallv(recv_buf_T,my_recv_count,my_recv_disp,MPI_T,
                  send_buf_T,my_send_count,my_send_disp,MPI_T,mpi_comm) ;


    // now unpack on the distributed side ...
    /*
    ii = 0;
    while ( ii < my_send_count_sum ) {

      int this_ida_global = send_buf_T[ii++] ; assert(this_ida_global < 0);
      this_ida_global     = -this_ida_global-1;

      int my_ida          = ixoxg[this_ida_global]; //XXX see note at top of page..

      int nod    = send_buf_T[ii++] ;

      // assumes that you already built my_i_i, so that
      // we can do some error checking ...
      if ( my_i_i[my_ida+1] - my_i_i[my_ida] != nod) {
        if ( mpi_rank == 0 ) {
          cout << my_ida << "    " << my_i_i[my_ida] << "   " << my_i_i[my_ida+1] << "   " << nod << endl ;
        }
      }

      assert(my_i_i[my_ida+1] - my_i_i[my_ida] == nod) ;

      for (int jj=0; jj < nod ; ++jj) {
        my_i_v[my_i_i[my_ida]+jj] = send_buf_T[ii++];
      }
    }
    */

    int* tmp_i = new int[send_count_sum+1];
    {
      ii = 0 ;
      int is = 0 ;

      while ( ii < my_send_count_sum ) {

        tmp_i[is] = ii ;
        T this_ida_global = send_buf_T[ii++];
        assert ( this_ida_global < 0);
        this_ida_global = -this_ida_global -1;

        T nod = send_buf_T[ii++];
        for (int jj = 0 ; jj < int(nod); ++jj)
          ++ii;

        is++;
      }
      assert ( is == send_count_sum) ;
    }

    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      ii           = tmp_i[is];

      T this_ida_global   = send_buf_T[ii++];
      assert( this_ida_global < 0) ;
      this_ida_global     = -this_ida_global-1;
      assert ( ix_global[ida] == this_ida_global) ;

      T nod    = send_buf_T[ii++];

      // assumes that you already built my_i_i, so that
      // we can do some error checking ...
      if ( my_i_i[ida+1] - my_i_i[ida] != nod) {
        if ( mpi_rank == 0 ) {
          cout << ida << "    " << my_i_i[ida] << "   " << my_i_i[ida+1] << "   " << nod << endl ;
        }
      }

      assert(my_i_i[ida+1] - my_i_i[ida] == nod) ;

      for (int jj=0; jj < nod ; ++jj) {
        my_i_v[my_i_i[ida]+jj] = send_buf_T[ii++];
      }
    }

    DELETE(tmp_i);
    DELETE(recv_buf_T);
    DELETE(send_buf_T);
    DELETE(my_recv_count);
    DELETE(my_recv_disp);
    DELETE(my_send_count);
    DELETE(my_send_disp);
  }

#endif

  // alot of trouble with this one...

  void pull_csr(const int * const my_i_i,uint8 * my_i_v,const int * const i_i, const uint8 * const i_v,const int8 * const ix_global,
                const int8 * const daora) {

    // XXX .. i need to build a map of my local indices to the global ones...
    // there may be a more efficient way of doing this ..
    /*
    std::map<int8,int> ixoxg ;
    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int8 this_ida_global = ix_global[ida];
      ixoxg[this_ida_global]     = ida ;
    }
    */

    MPI_Datatype MPI_T = MPI_UINT8; //MpiStuff::getMpiDatatype<T>();

    int * my_recv_count = new int[mpi_size];
    int * my_recv_disp  = new int[mpi_size];

    int * my_send_count = new int[mpi_size];
    int * my_send_disp  = new int[mpi_size];

    FOR_RANK my_recv_count[rank] = 0;

    int* recv_buf = new int[recv_count_sum];
    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida   = recv_index[ir];
      recv_buf[ir]    = int(i_i[ida+1]-i_i[ida])+2; // take some bookkeeping info too...
    }

    FOR_RECVNBR {

      const int offset   = recvIt->offset ;
      const int count    = recvIt->n ;
      const int nbr_rank = recvIt->nbr_rank ;

      for (int i =0; i < count ; ++i) {
        my_recv_count[nbr_rank] += recv_buf[i+offset];
      }
    }//recvnbr

    // build the send side ...
    MPI_Alltoall(my_recv_count,1,MPI_INT,my_send_count,1,MPI_INT,mpi_comm);

    // build the displacements ..
    int my_recv_count_sum = 0;
    int my_send_count_sum = 0;

    my_recv_disp[0] = 0;
    my_send_disp[0] = 0;
    for (int rank =1 ; rank < mpi_size ; ++rank) {
      my_recv_disp[rank] = my_recv_disp[rank-1] + my_recv_count[rank-1];
      my_recv_count_sum += my_recv_count[rank-1];

      my_send_disp[rank] = my_send_disp[rank-1] + my_send_count[rank-1];
      my_send_count_sum += my_send_count[rank-1];
    }

    my_recv_count_sum += my_recv_count[mpi_size-1];
    my_send_count_sum += my_send_count[mpi_size-1];

    DELETE(recv_buf);

    // pack ...
    uint8* recv_buf_T = new uint8[my_recv_count_sum];

    int ii = 0;
    for (int ir = 0 ; ir < recv_count_sum ; ++ir) {
      const int ida = recv_index[ir];

      // global index, the count, then the ranks ..
      int8 this_ida_global   = daora[mpi_rank] + ida;
      recv_buf_T[ii++] = this_ida_global;
      recv_buf_T[ii++] = (i_i[ida+1]-i_i[ida]) ;

      for (int dod = i_i[ida]; dod != i_i[ida+1]; ++dod) {
        recv_buf_T[ii++] = i_v[dod];
      }
    }//ir

    assert ( ii == my_recv_count_sum ) ;

    uint8* send_buf_T = new uint8[my_send_count_sum];

    // do the exchange
    MPI_Alltoallv(recv_buf_T,my_recv_count,my_recv_disp,MPI_T,
                  send_buf_T,my_send_count,my_send_disp,MPI_T,mpi_comm) ;


    // now unpack on the distributed side ...
    /*
    ii = 0;
    while ( ii < my_send_count_sum ) {

      int this_ida_global = send_buf_T[ii++] ; assert(this_ida_global < 0);
      this_ida_global     = -this_ida_global-1;

      int my_ida          = ixoxg[this_ida_global]; //XXX see note at top of page..

      int nod    = send_buf_T[ii++] ;

      // assumes that you already built my_i_i, so that
      // we can do some error checking ...
      if ( my_i_i[my_ida+1] - my_i_i[my_ida] != nod) {
        if ( mpi_rank == 0 ) {
          cout << my_ida << "    " << my_i_i[my_ida] << "   " << my_i_i[my_ida+1] << "   " << nod << endl ;
        }
      }

      assert(my_i_i[my_ida+1] - my_i_i[my_ida] == nod) ;

      for (int jj=0; jj < nod ; ++jj) {
        my_i_v[my_i_i[my_ida]+jj] = send_buf_T[ii++];
      }
    }
    */

    int* tmp_i = new int[send_count_sum+1];
    {
      ii = 0 ;
      int is = 0 ;
      while ( ii < my_send_count_sum ) {
	tmp_i[is++] = ii ;
        //uint8 this_ida_global = send_buf_T[ii++];
	++ii;
	const int nod = send_buf_T[ii++];
        ii += nod;
      }
      assert ( is == send_count_sum) ;
    }

    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int is = send_index[ida];
      ii           = tmp_i[is];

      int8 this_ida_global   = send_buf_T[ii++];
      assert ( ix_global[ida] == this_ida_global) ;

      const int nod    = send_buf_T[ii++];

      // assumes that you already built my_i_i, so that
      // we can do some error checking ...
      if ( my_i_i[ida+1] - my_i_i[ida] != nod) {
        if ( mpi_rank == 0 ) {
          cout << ida << "    " << my_i_i[ida] << "   " << my_i_i[ida+1] << "   " << nod << endl ;
	  assert(0);
        }
      }

      assert(my_i_i[ida+1] - my_i_i[ida] == nod) ;

      for (int jj=0; jj < nod ; ++jj) {
        my_i_v[my_i_i[ida]+jj] = send_buf_T[ii++];
      }
    }

    DELETE(tmp_i);
    DELETE(recv_buf_T);
    DELETE(send_buf_T);
    DELETE(my_recv_count);
    DELETE(my_recv_disp);
    DELETE(my_send_count);
    DELETE(my_send_disp);

  }

  void pull(int * my_i0,const int * i0) {

    int * recv_buf = new int[recv_count_sum];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = i0[ida];
    }

    int * send_buf = new int[send_count_sum];
    doReverseExchange<int>(recv_buf,send_buf);
    delete[] recv_buf;

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      my_i0[ida] = send_buf[is];
    }

    delete[] send_buf;

  }
 
  void pull(uint2 * my_i0,const uint2 * i0) {
  
    uint2 * recv_buf = new uint2[recv_count_sum];
    
    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = i0[ida];
    }
    
    uint2 * send_buf = new uint2[send_count_sum];
    doReverseExchange<uint2>(recv_buf,send_buf); 
    delete[] recv_buf;

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      my_i0[ida] = send_buf[is];
    }
    
    delete[] send_buf;
    
  }

  void pull(int8 * my_i0,const int8 * i0) {

    int8 * recv_buf = new int8[recv_count_sum];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = i0[ida];
    }

    int8 * send_buf = new int8[send_count_sum];
    doReverseExchange<int8>(recv_buf,send_buf) ;
    delete[] recv_buf;

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      my_i0[ida] = send_buf[is];
    }

    delete[] send_buf;
  }

  void pull(uint8 * my_i0,const uint8 * i0) {

    uint8 * recv_buf = new uint8[recv_count_sum];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = i0[ida];
    }

    uint8 * send_buf = new uint8[send_count_sum];
    doReverseExchange<uint8>(recv_buf,send_buf) ;
    delete[] recv_buf;

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      my_i0[ida] = send_buf[is];
    }

    delete[] send_buf;
  }

  void pull(int (*my_i0)[2],const int (*i0)[2]) {

    int * recv_buf = new int[recv_count_sum*2];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I2 recv_buf[2*ir+i] = i0[ida][i];
    }

    FOR_SENDNBR {
      sendIt->n *= 2 ;
      sendIt->offset *= 2 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 2 ;
      recvIt->offset *= 2 ;
    }


    int * send_buf = new int[send_count_sum*2];
    doReverseExchange<int>(recv_buf,send_buf) ;
    delete[] recv_buf;

    FOR_SENDNBR {
      sendIt->n /= 2 ;
      sendIt->offset /= 2 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 2 ;
      recvIt->offset /= 2 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I2 my_i0[ida][i] = send_buf[2*is+i];
    }

    delete[] send_buf;
  }


  void pull(int8 (*my_i0)[2],const int8 (*i0)[2]) {

    int8 * recv_buf = new int8[recv_count_sum*2];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I2 recv_buf[2*ir+i] = i0[ida][i];
    }
    
    FOR_SENDNBR {
      sendIt->n *= 2 ;
      sendIt->offset *= 2 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 2 ;
      recvIt->offset *= 2 ;
    }

    int8 * send_buf = new int8[send_count_sum*2];
    doReverseExchange<int8>(recv_buf,send_buf) ;
    delete[] recv_buf;

    FOR_SENDNBR {
      sendIt->n /= 2 ;
      sendIt->offset /= 2 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 2 ;
      recvIt->offset /= 2 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I2 my_i0[ida][i] = send_buf[2*is+i];
    }

    delete[] send_buf;

  }

  void pull(int8 (*my_i0)[3],const int8 (*i0)[3]) {

    int8 * recv_buf = new int8[recv_count_sum*3];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I3 recv_buf[3*ir+i] = i0[ida][i];
    }

    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }


    int8 * send_buf = new int8[send_count_sum*3];
    doReverseExchange<int8>(recv_buf,send_buf) ;
    delete[] recv_buf;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I3 my_i0[ida][i] = send_buf[3*is+i];
    }

    delete[] send_buf;

  }

  void pull(double * my_v0,const double * v0) {

    double * recv_buf = new double[recv_count_sum];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = v0[ida];
    }

    double * send_buf = new double[send_count_sum];
    doReverseExchange<double>(recv_buf,send_buf) ;
    delete[] recv_buf;

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      my_v0[ida] = send_buf[is];
    }

    delete[] send_buf;
  }

  void pull(double (*my_v0)[3],const double (*v0)[3]) {

    double * recv_buf = new double[recv_count_sum*3];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I3 recv_buf[3*ir+i] = v0[ida][i];
    }

    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }

    double * send_buf = new double[send_count_sum*3];
    doReverseExchange<double>(recv_buf,send_buf) ;
    delete[] recv_buf ;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I3 my_v0[ida][i] = send_buf[3*is+i];
    }

    delete[] send_buf;

  }

  void pull(double (*my_v0)[3][3],const double (*v0)[3][3]) {
    
    double * recv_buf = new double[recv_count_sum*9];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I3 FOR_J3 recv_buf[9*ir+3*i+j] = v0[ida][i][j];
    }

    FOR_SENDNBR {
      sendIt->n *= 9 ;
      sendIt->offset *= 9 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 9 ;
      recvIt->offset *= 9 ;
    }

    double * send_buf = new double[send_count_sum*9];
    doReverseExchange<double>(recv_buf,send_buf) ;
    delete[] recv_buf ;

    FOR_SENDNBR {
      sendIt->n /= 9 ;
      sendIt->offset /= 9 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 9 ;
      recvIt->offset /= 9 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I3 FOR_J3 my_v0[ida][i][j] = send_buf[9*is+3*i+j];
    }

    delete[] send_buf;

  }

  void pull(int (*my_v0)[3],const int (*v0)[3]) {

    int * recv_buf = new int[recv_count_sum*3];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I3 recv_buf[3*ir+i] = v0[ida][i];
    }

    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }

    int * send_buf = new int[send_count_sum*3];
    doReverseExchange<int>(recv_buf,send_buf) ;
    delete[] recv_buf ;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }

    for (int ida = 0; ida < my_nda; ++ida) {
      const int is = send_index[ida];
      FOR_I3 my_v0[ida][i] = send_buf[3*is+i];
    }

    delete[] send_buf;

  }

protected:

  template<class T>
  void doExchange(T * send_buf, T * recv_buf) {

#ifdef DDE_ALLTOALL

    //if ( mpi_rank == 0 )
    //  cout << " > Using a2a dde " << endl ;


    int * send_count = new int[mpi_size];
    int * recv_count = new int[mpi_size];
    int * send_disp  = new int[mpi_size];
    int * recv_disp  = new int[mpi_size];

    FOR_RANK {
       send_count[rank]  = 0;
       recv_count[rank]  = 0;
    }

    FOR_SENDNBR send_count[sendIt->nbr_rank] = sendIt->n ;
    FOR_RECVNBR recv_count[recvIt->nbr_rank] = recvIt->n ;

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>();

    //build the displacements from the offset values..
    FOR_RANK {
       send_disp[rank] = -1;
       recv_disp[rank] = -1;
    }

    send_disp[0] = 0;
    recv_disp[0] = 0;
    FOR_SENDNBR {
      send_disp[sendIt->nbr_rank] = sendIt->offset;
    }

    FOR_RECVNBR {
      recv_disp[recvIt->nbr_rank] = recvIt->offset;
    }

    assert ( recv_disp[0] == 0 ) ;
    assert ( send_disp[0] == 0 ) ;

    for (int i=1 ; i < mpi_size ; ++i) {
      if ( recv_disp[i] == -1) {
        assert ( recv_count[i] == 0 );
        recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
      }

      if ( send_disp[i] == -1) {
        assert ( send_count[i]  == 0 );
        send_disp[i] = send_count[i-1] + send_disp[i-1];
      }
    }//i


    MPI_Alltoallv( send_buf, send_count, send_disp, MPI_T,
                   recv_buf, recv_count, recv_disp, MPI_T, mpi_comm);

    delete[] recv_disp ;
    delete[] send_disp ;
    delete[] send_count ;
    delete[] recv_count ;

#else

    const int nreq_max = recvNbrVec.size(); //sendNbrVec.size();
    //const int nreq_max = sendNbrVec.size();
    MPI_Request * reqArr = new MPI_Request[nreq_max];
    MPI_Status  * statArr = new MPI_Status[nreq_max];
    int req = 0 ;

    //vector<PeerNbr>::iterator mySend ;
    //vector<PeerNbr>::iterator myRecv ;
    PeerNbr* mySend  ;
    PeerNbr* myRecv  ;

    bool isSelfSend = false ;


   // if ( mpi_rank == 0 )
    //  cout << " > Using pt2pt dde " << endl ;

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>() ;

    /*
    FOR_SENDNBR {

      int offset   = sendIt->offset ;
      int count    = sendIt->n ;
      int nbr_rank = sendIt->nbr_rank ;

      if ( nbr_rank != mpi_rank ) {

        MPI_Issend(&send_buf[offset], count, MPI_T, nbr_rank, mpi_rank, mpi_comm,
                   &reqArr[req]) ;
        ++req;
      }
      else {
        isSelfSend = true ;
        mySend = sendIt ;
      }
    }


    // block on the recv end ..
    FOR_RECVNBR {

      int offset   = recvIt->offset ;
      int count    = recvIt->n ;
      int nbr_rank = recvIt->nbr_rank ;

      MPI_Status status ;
      if ( nbr_rank != mpi_rank ) {
        MPI_Recv(&recv_buf[offset], count, MPI_T, nbr_rank, nbr_rank, mpi_comm,&status) ;
      }
      else {
        assert(isSelfSend) ;
        myRecv = recvIt ;
      }
    }
    */

    MPI_Request* sendReq = new MPI_Request[sendNbrVec.size()];
    int sreq = 0;

    FOR_RECVNBR {
      int offset   = recvIt->offset ;
      int count    = recvIt->n ;
      int nbr_rank = recvIt->nbr_rank ;

      MPI_Status status ;
      if ( nbr_rank != mpi_rank ) {
        MPI_Irecv(&recv_buf[offset], count, MPI_T, nbr_rank, nbr_rank, mpi_comm,&reqArr[req]) ;
        ++req;
      }
      else {
        isSelfSend = true ;
        myRecv = recvIt ;
      }
    }

    FOR_SENDNBR {

      int offset   = sendIt->offset ;
      int count    = sendIt->n ;
      int nbr_rank = sendIt->nbr_rank ;

      if ( nbr_rank != mpi_rank ) {
	MPI_Ssend(&send_buf[offset], count, MPI_T, nbr_rank, mpi_rank, mpi_comm);
	//MPI_Issend(&send_buf[offset],count,MPI_T,nbr_rank,mpi_rank,mpi_comm,&sendReq[sreq]);
        ++sreq ;
      }
      else {
        assert( isSelfSend) ;
        mySend = sendIt ;
      }
    }




    if ( isSelfSend ) {
      int send_offset = mySend->offset ;
      int send_count  = mySend->n ;
      int recv_offset = myRecv->offset ;
      int recv_count  = myRecv->n ;

      assert(send_count == recv_count) ;

      for (int i = 0 ; i < send_count ; ++i) {
        recv_buf[i+recv_offset] = send_buf[i+send_offset];
      }
    }

    // check to make sure that everyone finished, but you shouldnt have to...
    if ( req > 0 ) {
      MPI_Waitall(req, reqArr, statArr) ;
    }

    delete[] reqArr ;
    delete[] statArr;
    delete[] sendReq;

#endif

  }


  // recv--->send (recv,send are based on the way that the peerNbrs were built)
  template<class T>
  void doReverseExchange(T * recv_buf, T * send_buf) {

#ifdef DDE_ALLTOALL

    //if ( mpi_rank == 0 )
    //  cout << " > Using a2a dde " << endl ;

    int * send_count = new int[mpi_size];
    int * recv_count = new int[mpi_size];
    int * send_disp  = new int[mpi_size];
    int * recv_disp  = new int[mpi_size];

    FOR_RANK {
       send_count[rank]  = 0;
       recv_count[rank]  = 0;
    }

    FOR_SENDNBR send_count[sendIt->nbr_rank] = sendIt->n ;
    FOR_RECVNBR recv_count[recvIt->nbr_rank] = recvIt->n ;

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>();

    //build the displacements from the offset values..
    FOR_RANK {
       send_disp[rank] = -1;
       recv_disp[rank] = -1;
    }

    send_disp[0] = 0;
    recv_disp[0] = 0;
    FOR_SENDNBR {
      send_disp[sendIt->nbr_rank] = sendIt->offset;
    }

    FOR_RECVNBR {
      recv_disp[recvIt->nbr_rank] = recvIt->offset;
    }

    assert ( recv_disp[0] == 0 ) ;
    assert ( send_disp[0] == 0 ) ;

    for (int i=1 ; i < mpi_size ; ++i) {
      if ( recv_disp[i] == -1) {
        assert ( recv_count[i] == 0 );
        recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
      }

      if ( send_disp[i] == -1) {
        assert ( send_count[i]  == 0 );
        send_disp[i] = send_count[i-1] + send_disp[i-1];
      }
    }//i


    MPI_Alltoallv( recv_buf, recv_count, recv_disp, MPI_T,
                   send_buf, send_count, send_disp, MPI_T, mpi_comm);

    delete[] recv_disp ;
    delete[] send_disp ;
    delete[] send_count ;
    delete[] recv_count ;

#else

    PeerNbr* mySend ;
    PeerNbr* myRecv ;

    //if ( mpi_rank == 0 )
    //  cout << " > Using pt2pt dde " << endl ;

    MPI_Datatype MPI_T = MpiStuff::getMpiDatatype<T>() ;

    bool isSelfSend = false ;
    const int nreq_max   = sendNbrVec.size(); //recvNbrVec.size() ;
    MPI_Request * reqArr = new MPI_Request[nreq_max];
    MPI_Status * statArr = new MPI_Status[nreq_max];
    int req = 0 ;

    /*
    FOR_RECVNBR {

      int offset   = recvIt->offset ;
      int count    = recvIt->n ;
      int nbr_rank = recvIt->nbr_rank ;

      if ( nbr_rank != mpi_rank ) {
        MPI_Issend(&recv_buf[offset], count, MPI_T, nbr_rank, mpi_rank, mpi_comm,
                   &reqArr[req]) ;
        ++req;
      }
      else {
        isSelfSend = true ;
        myRecv = recvIt ;
      }
    }


    // block on the recv end ..
    FOR_SENDNBR {

      int offset   = sendIt->offset ;
      int count    = sendIt->n ;
      int nbr_rank = sendIt->nbr_rank ;

      MPI_Status status ;
      if ( nbr_rank != mpi_rank ) {
        MPI_Recv(&send_buf[offset], count, MPI_T, nbr_rank, nbr_rank, mpi_comm, &status) ;
      }
      else {
        assert(isSelfSend);
        mySend = sendIt;
      }
    }
    */

    MPI_Request* recvReq = new MPI_Request[recvNbrVec.size()];
    int rreq = 0;

    FOR_SENDNBR {

      int offset   = sendIt->offset ;
      int count    = sendIt->n ;
      int nbr_rank = sendIt->nbr_rank ;

      MPI_Status status ;
      if ( nbr_rank != mpi_rank ) {
        MPI_Irecv(&send_buf[offset], count, MPI_T, nbr_rank, nbr_rank, mpi_comm, &reqArr[req]) ;
        ++req;
      }
      else {
        isSelfSend = true;
        mySend = sendIt;
      }
    }

    FOR_RECVNBR {

      int offset   = recvIt->offset ;
      int count    = recvIt->n ;
      int nbr_rank = recvIt->nbr_rank ;

      if ( nbr_rank != mpi_rank ) {
        MPI_Ssend(&recv_buf[offset], count, MPI_T, nbr_rank, mpi_rank, mpi_comm);
        //MPI_Issend(&recv_buf[offset],count,MPI_T,nbr_rank,mpi_rank,mpi_comm,&recvReq[rreq]);
        ++rreq;
      }
      else {
        assert(isSelfSend) ;
        myRecv = recvIt ;
      }
    }


    if ( isSelfSend) {
      int send_offset = mySend->offset ;
      int send_count  = mySend->n ;
      int recv_offset = myRecv->offset ;
      int recv_count  = myRecv->n ;

      assert(send_count == recv_count) ;

      for (int i = 0 ; i < send_count ; ++i) {
        send_buf[i+send_offset] = recv_buf[i+recv_offset];
      }
    }


    // check to make sure that everyone finished, but you shouldnt have to...
    if ( req > 0 )
      MPI_Waitall(req, reqArr, statArr) ;

    delete[] reqArr ;
    delete[] statArr ;
    delete[] recvReq;
#endif
  }


  void init(const int * my_ida_global,const int my_nda,const int daora[]) {


    //COUT1(" > Starting dd exchanger init (int)... " ) ;

    send_count = NULL ;
    send_disp  = NULL ;
    send_index = NULL ;

    recv_count = NULL ;
    recv_disp  = NULL ;
    recv_index = NULL ;

    int * recv_count = new int[mpi_size];
    FOR_RANK recv_count[rank] = 0 ;

    // here we have some data that knows where it needs to go in terms of some distributed
    // arrays across the processors...

    this->my_nda = my_nda;
    vector<IntPair> intPairVec(my_nda);
    for (int ida = 0; ida < my_nda; ++ida) {
      // check that range of data index is valid...
      assert((my_ida_global[ida] >= daora[0])&&(my_ida_global[ida] < daora[mpi_size]));
      // populate the int pair...
      intPairVec[ida].i1 = my_ida_global[ida];
      intPairVec[ida].i2 = ida;
    }

    // and stable sort based on my_ida_global, then ida. Here we prefer a stable sort
    // because it simplifies certain actions when packing the send side...

    std::sort(intPairVec.begin(),intPairVec.end(),IntPairCompare12());

    // and figure out where to send the data...
    sendNbrVec.clear() ;
    int rank = -1;
    int current_ida_global = -1;
    send_count_sum = 0 ;
    PeerNbr* prNbr = NULL ;

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;

    int * send_buf_int = new int[my_nda] ;

    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int ida_global = intPairVec[ida].i1 ;
      if ( ida_global != current_ida_global ) {
        assert ( ida_global > current_ida_global) ;
        if ( ida_global >= daora[rank+1] ) {

          while ( ida_global >= daora[rank+1])
            ++rank ;

          sendNbrVec.push_front(PeerNbr(rank,0,send_count_sum)) ;
          prNbr   = sendNbrVec.front() ;

        }


        // find the send count of the last entry in the vector ...
        assert( prNbr != NULL ) ;
        prNbr->n++;
        send_count[rank] += 1 ;
        send_buf_int[send_count_sum++] = ida_global ;
        current_ida_global = ida_global ;
      }
      intPairVec[ida].i1 = send_count_sum-1;
    }//ida ...

    assert(int(intPairVec.size()) == my_nda);
    send_index = new int[my_nda];
    for (int ida = 0; ida < my_nda; ++ida)
      send_index[intPairVec[ida].i2] = intPairVec[ida].i1;

    intPairVec.clear();

    //MPI_Win_fence(0,recvWin) ;
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    DELETE(send_count) ;


    recvNbrVec.clear() ;
    recv_count_sum = 0 ;
    FOR_RANK {
      if ( recv_count[rank] > 0 ) {
        //recvNbrVec.push_back(PeerNbr(rank,recv_count[rank],recv_count_sum)) ;
        recvNbrVec.push_front(PeerNbr(rank,recv_count[rank], recv_count_sum)) ;
        recv_count_sum += recv_count[rank];
      }
    }//rank ...


    //DELETE(send_count) ;
    DELETE(recv_count) ;


    int * recv_buf_int = new int[recv_count_sum] ;
    doExchange<int>(send_buf_int, recv_buf_int ) ;


    // set nda...

    nda = daora[mpi_rank+1] - daora[mpi_rank];

    // build the recv_index. Note that da_flag is for checking here for the case when
    // all the ditributed data should be set once and only once...

    // here da_flag is just for checking and can be discarded...

    /*
    int * da_flag = new int[nda];
    for (int ida = 0; ida < nda; ++ida)
      da_flag[ida] = 0;
    */

    recv_index = new int[recv_count_sum];
    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int this_ida_global = recv_buf_int[ir];
      assert((this_ida_global >= daora[mpi_rank])&&(this_ida_global < daora[mpi_rank+1]));
      const int this_ida = this_ida_global - daora[mpi_rank];
      recv_index[ir] = this_ida;
      //da_flag[this_ida] += 1;
    }

    DELETE(send_buf_int);
    DELETE(recv_buf_int);

    // it is not neccessary, but we can check that each member in the data was
    // touched once...

    /*
    for (int ida = 0; ida < nda; ++ida)
      assert(da_flag[ida] == 1);
    delete[] da_flag;
    */

    //COUT1(" > Finishing dd exchanger init (int)... " ) ;
  }

  void init(const int8 * my_ida_global,const int my_nda,const int8 daora[]) {

    //COUT1(" > Starting dd exchanger init (int8)... " ) ;

    send_count = NULL ;
    send_disp  = NULL ;
    send_index = NULL ;

    recv_count = NULL ;
    recv_disp  = NULL ;
    recv_index = NULL ;

    int * recv_count = new int[mpi_size];

    FOR_RANK recv_count[rank] = 0 ;

    // here we have some data that knows where it needs to go in terms of some distributed
    // arrays across the processors...

    this->my_nda = my_nda;
    vector<Int8IntPair> intPairVec(my_nda);
    for (int ida = 0; ida < my_nda; ++ida) {
      // check that range of data index is valid...
      //cout << my_ida_global[ida] << " " << daora[0] << " " << daora[mpi_size] << endl;
      assert((my_ida_global[ida] >= daora[0])&&(my_ida_global[ida] < daora[mpi_size]));
      // populate the int pair...
      intPairVec[ida].i1 = my_ida_global[ida];
      intPairVec[ida].i2 = ida;
    }

    // and stable sort based on my_ida_global, then ida. Here we prefer a stable sort
    // because it simplifies certain actions when packing the send side...

    std::sort(intPairVec.begin(),intPairVec.end(),Int8IntPairCompare12());

    // and figure out where to send the data...
    sendNbrVec.clear() ;
    int rank = -1;
    int8 current_ida_global = -1;
    send_count_sum = 0 ;
    PeerNbr* prNbr = NULL ;

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0 ;

    int8 * send_buf_int = new int8[my_nda] ;

    for (int ida = 0 ; ida < my_nda ; ++ida) {
      const int8 ida_global = intPairVec[ida].i1 ;
      if ( ida_global != current_ida_global ) {
        assert ( ida_global > current_ida_global) ;
        if ( ida_global >= daora[rank+1] ) {

          while ( ida_global >= daora[rank+1])
            ++rank ;

          sendNbrVec.push_front(PeerNbr(rank,0,send_count_sum)) ;
          prNbr = sendNbrVec.front();

        }

        // find the send count of the last entry in the vector ...
        assert( prNbr != NULL ) ;
        send_count[rank] += 1;
        prNbr->n++;
        send_buf_int[send_count_sum++] = ida_global ;
        current_ida_global = ida_global ;
      }
      intPairVec[ida].i1 = int(send_count_sum-1);
    }//ida ...

    assert(int(intPairVec.size()) == my_nda);
    send_index = new int[my_nda];
    for (int ida = 0; ida < my_nda; ++ida)
      send_index[intPairVec[ida].i2] = int(intPairVec[ida].i1);

    intPairVec.clear();

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    DELETE(send_count) ;

    recvNbrVec.clear() ;
    recv_count_sum = 0 ;
    FOR_RANK {
      if ( recv_count[rank] > 0 ) {
        //recvNbrVec.push_back(PeerNbr(rank,recv_count[rank],recv_count_sum)) ;
        recvNbrVec.push_front(PeerNbr(rank,recv_count[rank],recv_count_sum)) ;
        recv_count_sum += recv_count[rank];
      }
    }//rank ...


    //DELETE(send_count) ;
    //DELETE(recv_count) ;


    int8 * recv_buf_int = new int8[recv_count_sum] ;
    doExchange<int8>(send_buf_int, recv_buf_int ) ;

    DELETE(send_buf_int);

    // set nda...

    nda = daora[mpi_rank+1] - daora[mpi_rank];

    // build the recv_index. Note that da_flag is for checking here for the case when
    // all the ditributed data should be set once and only once...

    // here da_flag is just for checking and can be discarded...

    /*
    int * da_flag = new int[nda];
    for (int ida = 0; ida < nda; ++ida)
      da_flag[ida] = 0;
    */

    recv_index = new int[recv_count_sum];
    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int8 this_ida_global = recv_buf_int[ir];
      assert((this_ida_global >= daora[mpi_rank])&&(this_ida_global < daora[mpi_rank+1]));
      const int this_ida = int(this_ida_global - daora[mpi_rank]);
      recv_index[ir] = this_ida;
      //da_flag[this_ida] += 1;
    }

    DELETE(recv_buf_int) ;

    // it is not neccessary, but we can check that each member in the data was
    // touched once...

    /*
    for (int ida = 0; ida < nda; ++ida)
      assert(da_flag[ida] == 1);
    delete[] da_flag;
    */

    DELETE(recv_count) ;

    //COUT1(" > Finishing dd exchanger init (int8)... " ) ;
  }

  void clearSendIndexIfPossible() {

    // in some cases, the send buffer packs the same order as the
    // passed data buffer. In the future, we can detect this and
    // modify the front-end operations accordingly...

    /*
    int my_clear = 1;
    for (int ida = 0; ida < my_nda; ++ida) {
      if (send_index[ida] != ida) {
	my_clear = 0;
	break;
      }
    }

    int clear;
    MPI_Reduce(&my_clear,&clear,1,MPI_INT,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > DistributedDataExchanger: could clear " << clear << " out of " << mpi_size << " send_buffers." << endl;
    */

  }

};

// ========================================================================
// A DistributedDataExchangerReducer builds on the class above to reduce
// following the pull operation...
// ========================================================================

class DistributedDataExchangerReducer : public DistributedDataExchanger {

public:

  DistributedDataExchangerReducer(int * my_ida_global,const int my_nda,const int daora[]) :
    DistributedDataExchanger(my_ida_global,my_nda,daora) {

    // now convert the my_ida_global to the unpack index in
    // is should be the same array used as daoda2_v below...

    for (int ida = 0; ida < my_nda; ++ida) {
      my_ida_global[ida] = send_index[ida];
    }

    DELETE(send_index);
  }

  void pullAndReduce(double * da2,const int nda2,const int * daoda2_i,const int * daoda2_v,const double * daoda2_wgt,
		     const double * v0) {

    double * recv_buf = new double[recv_count_sum];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      recv_buf[ir] = v0[ida];
    }

    double * send_buf = new double[send_count_sum];
    doReverseExchange<double>(recv_buf,send_buf);

    for (int ida2 = 0; ida2 < nda2; ++ida2) {
      da2[ida2] = 0.0;
      for (int dod = daoda2_i[ida2]; dod != daoda2_i[ida2+1]; ++dod) {
	const int ida = daoda2_v[dod];
	da2[ida2] += daoda2_wgt[dod]*send_buf[ida];
      }
    }

    delete[] send_buf;
  }

  void pullAndReduce(double (*da2)[3],const int nda2,const int * daoda2_i,const int * daoda2_v,const double * daoda2_wgt,
		     const double (*v0)[3]) {


    double * recv_buf = new double[recv_count_sum*3];

    for (int ir = 0; ir < recv_count_sum; ++ir) {
      const int ida = recv_index[ir];
      FOR_I3 recv_buf[3*ir+i] = v0[ida][i];
    }


    FOR_SENDNBR {
      sendIt->n *= 3 ;
      sendIt->offset *= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n *= 3 ;
      recvIt->offset *= 3 ;
    }

    double * send_buf = new double[send_count_sum*3];
    doReverseExchange<double>(recv_buf,send_buf);
    delete[] recv_buf;

    FOR_SENDNBR {
      sendIt->n /= 3 ;
      sendIt->offset /= 3 ;
    }

    FOR_RECVNBR {
      recvIt->n /= 3 ;
      recvIt->offset /= 3 ;
    }

    for (int ida2 = 0; ida2 < nda2; ++ida2) {
      FOR_I3 da2[ida2][i] = 0.0;
      for (int dod = daoda2_i[ida2]; dod != daoda2_i[ida2+1]; ++dod) {
	const int ida = daoda2_v[dod];
	FOR_I3 da2[ida2][i] += daoda2_wgt[dod]*send_buf[3*ida+i];
      }
    }

    delete[] send_buf;
  }

  void pull() {
    CERR("a simple pull is not available in a DistributedDataExchangerReducer");
  }
};
#endif
