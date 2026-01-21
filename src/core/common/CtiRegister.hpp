#ifndef _CTI_REGISTER_HPP_
#define _CTI_REGISTER_HPP_

#include "Common.hpp"
#include "CommonIo.hpp"
#include "MpiStuff.hpp"
using namespace MpiStuff;
using namespace std;

#include "DistributedDataExchanger.hpp"

// registration bits...
// there is definately a more efficient way to do this

// data is differentiated by type and topology. Type
// means value, scalar, vector, etc. Topology
// refers to the location of the data, like cv,
// bf. For the case of bf, where multiple bf's exist,
// there is an additional index specified in the most
// significant part of the bits. Topologies must be the
// same to operate on data.

// NOTE: if you add more here, then the index shift used in setTopology(topo,index) needs to be offset more
#define INDEX_SHIFT 20 // 2^(32-INDEX_SHIFT)-1 = 4095
#define TOPOLOGY_MASK ((4095<<INDEX_SHIFT)|CV_DATA|BF_DATA|NO_DATA|FA_DATA|EF_DATA|LP_DATA)

// read and/or writable
#define NO_READWRITE_DATA 0
#define READ_DATA         1
#define WRITE_DATA        2
#define READWRITE_DATA    3 // (READ_DATA|WRITE_DATA)
#define CAN_WRITE_DATA    4 // has a dde pointer
// topology (for IN*/DN*)
#define CV_DATA           8
#define BF_DATA           16 // zone is in bits>>INDEX_SHIFT
#define FA_DATA           32 // bool for signed is in bits>>INDEX_SHIFT (0 unsigned, 1 signed)
#define UNSIGNED_FA_DATA  FA_DATA
#define SIGNED_FA_DATA    1048608 // FA_DATA|(1<<INDEX_SHIFT)
#define EF_DATA           64
#define NO_DATA           128
#define LP_DATA           256
// data type
#define I_DATA            512    // integer (I0)
#define D_DATA            1024   // double (D0)
#define D3_DATA           2048   // double[3]
#define DN_DATA           4096   // double[N] (D1)
#define DN3_DATA          8192   // double[N][3] (D2)
#define DN33_DATA         16384  // double[N][3][3] (D3?)
#define IN_DATA           32768  // integer[N] (I1)
#define X_DATA            65536  // used to record x-data associated with a particular "n"...
#define L_DATA            131072 // used to record lengthscale data for a particular "n"...
#define S_DATA            262144 // string (S0)

// ==================================================================
// this is what it should look like...
// topology (for IN*/DN*)
//#define CV_DATA           1
//#define BF_DATA           2 // zone is in bits>>INDEX_SHIFT
//#define FA_DATA_           3 // bool for signed is in bits>>INDEX_SHIFT (0 unsigned, 1 signed)
//#define UNSIGNED_FA_DATA_  3
//#define SIGNED_FA_DATA_    4
//#define EF_DATA           5
//#define NO_DATA           6
//#define LP_DATA           7
// other topologies to 15...

// data types
//#define I_DATA            16 // integer (I0)
//#define D_DATA            17 // double (D0)
//#define D3_DATA           18 // double[3]
//#define DN_DATA           19 // double[N] (D1)
//#define DN3_DATA          20 // double[N][3] (D2)
//#define DN33_DATA         21 // double[N][3][3] (D3?)
//#define IN_DATA           22 // integer[N] (I1)
// other data types to 31...

// the remaining are bits...
//#define READ_DATA         (1<<5)  // 64
//#define WRITE_DATA        (1<<6)  // 128
//#define READWRITE_DATA    (3<<5)  // (READ_DATA|WRITE_DATA)
//#define CAN_WRITE_DATA    (1<<7)  // has a dde pointer
//#define X_DATA            (1<<8)  // used to record x-data associated with a particular "n"...
//#define L_DATA            (1<<9)  // used to record lengthscale data for a particular "n"...

namespace CtiRegister {

  void setEvalVerbosity(const bool b);

  enum CtiDataError {
    CTI_DATA_OK,
    CTI_DATA_ARG_COUNT,
    CTI_DATA_NOT_FOUND,
    CTI_DATA_NOT_VALID
  };

  class CtiData {
  private:

    int * size_ptr; // pointer to the number of data entries
    void * data_ptr; // pointer or pointer-pointer to data (needs to be casted)
    bool mem_flag; // do I manage the data?
    bool access_flag; // used for telling which type of access (different from mem_flag sometimes)
    uint bits; // read/write bits
    int stride; // (INTEGER) stride for structs/classes
    int flag; // used to indicate operations on registered data
    int offset; // (INTEGER) offset into struct/class

  public:

    // making these public bc the solver needs to check if these have been set
    // to allow data writing of unregisteredctidata

    // TODO: get rid of this eventually and rely on the NData capability

    int8 ** xora_ptr; // pointer-pointer to xora for i/o
    DistributedDataExchanger ** dde_ptr; // pointer-pointer to dde for i/o

  public:

    CtiData();
    CtiData(const CtiData& data);

    ~CtiData();

    void clear();

    bool empty() const;

    void reset();

    void copy(const CtiData& data);

    // ==================================
    // accessor functions...
    // ==================================

    void new_s() {
      assert(bits == 0);
      setType(S_DATA);
      assert(data_ptr == NULL);
      data_ptr = new string;
      assert(!mem_flag);
      mem_flag = true;
      access_flag = true;
    }

    void new_i() {
      assert(bits == 0);
      setType(I_DATA);
      assert(data_ptr == NULL);
      data_ptr = new int;
      assert(!mem_flag);
      mem_flag = true;
      access_flag = true;
    }

    void new_d() {
      assert(bits == 0);
      setType(D_DATA);
      assert(data_ptr == NULL);
      data_ptr = new double;
      assert(!mem_flag);
      mem_flag = true;
      access_flag = true;
    }

    void new_d3() {
      assert(bits == 0);
      setType(D3_DATA);
      assert(data_ptr == NULL);
      data_ptr = new double[3];
      assert(!mem_flag);
      mem_flag = true;
      access_flag = true;
    }

    // include optional size_ptr...
    // sometimes we want to delay alloc (e.g., register/init stats)

    // TODO: this delayed allocation adds complexity here and is only
    // used for stats. Should be promoted to some structure
    // that manages memory similar to the way solvers manage memory,
    // but for now...

    void new_dn(const CtiData& data,const bool b_init = true) {
      //if (n_ptr != NULL) assert(n_ptr == data.size_ptr);
      assert(bits == 0);
      assert(data.getType() >= DN_DATA);
      setType(DN_DATA);
      setTopology(data.getTopology());
      assert(size_ptr == NULL);
      assert(data.size_ptr != NULL);
      size_ptr = data.size_ptr;
      assert(data_ptr == NULL);
      if (b_init) {
        data_ptr = new double[*size_ptr];
      }
      stride = 2;
      mem_flag = true;
      access_flag = true;
      xora_ptr = data.xora_ptr;
      dde_ptr = data.dde_ptr;
    }

    void new_dn(const CtiData& data1,const CtiData& data2) {
      assert(bits == 0);
      assert(data1.getType() >= DN_DATA);
      assert(data2.getType() >= DN_DATA);
      setType(DN_DATA);
      const int unindexed_topo1 = data1.getUnindexedTopology();
      const int unindexed_topo2 = data2.getUnindexedTopology();
      assert(unindexed_topo1 == unindexed_topo2);
      if (unindexed_topo1 != FA_DATA) {
        assert(data1.getTopology() == data2.getTopology());
        setTopology(data1.getTopology());
      }
      else {
        // CI: what is this about???
        const int bit1 = data1.getIndex();
        const int bit2 = data2.getIndex();
        assert((bit1 == 0)||(bit1 == 1));
        assert((bit2 == 0)||(bit2 == 1));
        setTopology(unindexed_topo1,bit1 ^ bit2);
      }
      assert(data1.size() == data2.size());
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = data1.size();
      size_ptr = data1.size_ptr;
      assert(data_ptr == NULL);
      data_ptr = new double[*size_ptr];
      stride = 2;
      mem_flag = true;
      access_flag = true;
      setDdeStuff(data1,data2);
    }

    void new_dn(const int _topology,int &_n) {
      assert(bits == 0);
      setType(DN_DATA);
      setTopology(_topology);
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = _n;
      size_ptr = &_n;
      assert(data_ptr == NULL);
      data_ptr = new double[*size_ptr];
      stride = 2;
      mem_flag = true;
      access_flag = true;
      xora_ptr = NULL;
      dde_ptr = NULL;
    }

    void alloc_dn() {
      //int* tmp = size_ptr;
      //size_ptr = new int;
      //*size_ptr = *tmp;
      assert(size_ptr != NULL);
      assert(data_ptr == NULL);
      data_ptr = new double[*size_ptr];
    }

    // include optional size_ptr...
    void new_dn3(const CtiData& data,const bool b_init = true) {
      assert(bits == 0);
      assert(data.getType() >= DN_DATA);
      setType(DN3_DATA);
      setTopology(data.getTopology());
      assert(size_ptr == NULL);
      assert(data.size_ptr != NULL);
      size_ptr = data.size_ptr;
      assert(data_ptr == NULL);
      if (b_init) {
        data_ptr = new double[(*size_ptr)*3];
      }
      stride = 6;
      mem_flag = true;
      access_flag = true;
      // copy whether NULL or not...
      xora_ptr = data.xora_ptr;
      dde_ptr = data.dde_ptr;
    }
    void alloc_dn3() {
      //int* tmp = size_ptr;
      //size_ptr = new int;
      //*size_ptr = *tmp;
      assert(size_ptr != NULL);
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*3];
    }

    void new_dn3(const CtiData& data1,const CtiData& data2) {
      assert(bits == 0);
      assert(data1.getType() >= DN_DATA);
      assert(data2.getType() >= DN_DATA);
      setType(DN3_DATA);
      const int unindexed_topo1 = data1.getUnindexedTopology();
      const int unindexed_topo2 = data2.getUnindexedTopology();
      assert(unindexed_topo1 == unindexed_topo2);
      if (unindexed_topo1 != FA_DATA) {
        assert(data1.getTopology() == data2.getTopology());
        setTopology(data1.getTopology());
      }
      else {
        const int bit1 = data1.getIndex();
        const int bit2 = data2.getIndex();
        assert((bit1 == 0)||(bit1 == 1));
        assert((bit2 == 0)||(bit2 == 1));
        setTopology(unindexed_topo1,bit1 ^ bit2);
      }
      assert(data1.size() == data2.size());
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = data1.size();
      size_ptr = data1.size_ptr;
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*3];
      stride = 6;
      mem_flag = true;
      access_flag = true;
      setDdeStuff(data1,data2);
    }

    void new_dn3(const int _topology,int &_n) {
      assert(bits == 0);
      setType(DN3_DATA);
      setTopology(_topology);
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = _n;
      size_ptr = &_n;
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*3];
      stride = 6;
      mem_flag = true;
      access_flag = true;
      xora_ptr = NULL;
      dde_ptr = NULL;
    }

    // TODO: get rid of this delayed allocation
    void new_dn33(const CtiData& data,const bool b_init = true) {
      //if (n_ptr != NULL) assert(data.size_ptr == n_ptr);
      assert(bits == 0);
      assert(data.getType() >= DN_DATA);
      setType(DN33_DATA);
      setTopology(data.getTopology());
      assert(size_ptr == NULL);
      assert(data.size_ptr != NULL);
      size_ptr = data.size_ptr;
      assert(data_ptr == NULL);
      if (b_init) {
        data_ptr = new double[(*size_ptr)*9];
      }
      stride = 18;
      mem_flag = true;
      access_flag = true;
      // copy whether NULL or not...
      xora_ptr = data.xora_ptr;
      dde_ptr = data.dde_ptr;
    }

    void new_dn33(const CtiData& data1,const CtiData& data2) {
      assert(bits == 0);
      assert(data1.getType() >= DN_DATA);
      assert(data2.getType() >= DN_DATA);
      setType(DN33_DATA);
      const int unindexed_topo1 = data1.getUnindexedTopology();
      const int unindexed_topo2 = data2.getUnindexedTopology();
      assert(unindexed_topo1 == unindexed_topo2);
      if (unindexed_topo1 != FA_DATA) {
        assert(data1.getTopology() == data2.getTopology());
        setTopology(data1.getTopology());
      }
      else {
        const int bit1 = data1.getIndex();
        const int bit2 = data2.getIndex();
        assert((bit1 == 0)||(bit1 == 1));
        assert((bit2 == 0)||(bit2 == 1));
        setTopology(unindexed_topo1,bit1 ^ bit2);
      }
      assert(data1.size() == data2.size());
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = data1.size();
      size_ptr = data1.size_ptr;
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*9];
      stride = 18;
      mem_flag = true;
      access_flag = true;
      setDdeStuff(data1,data2);
    }

    void new_dn33(const int _topology,int &_n) {
      assert(bits == 0);
      setType(DN33_DATA);
      setTopology(_topology);
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = _n;
      size_ptr = &_n;
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*9];
      stride = 18;
      mem_flag = true;
      access_flag = true;
      xora_ptr = NULL;
      dde_ptr = NULL;
    }

    void alloc_dn33() {
      //int* tmp = size_ptr;
      //size_ptr = new int;
      //*size_ptr = *tmp;
      assert(size_ptr != NULL);
      assert(data_ptr == NULL);
      data_ptr = new double[(*size_ptr)*9];
    }

    void new_in(const CtiData& data) {
      assert(bits == 0);
      assert(data.getType() >= DN_DATA);
      setType(IN_DATA);
      setTopology(data.getTopology());
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = data.size();
      size_ptr = data.size_ptr;
      assert(data_ptr == NULL);
      data_ptr = new int[*size_ptr];
      stride = 1;
      mem_flag = true;
      access_flag = true;
      xora_ptr = data.xora_ptr;
      dde_ptr = data.dde_ptr;
    }

    void new_in(const CtiData& data1,const CtiData& data2) {
      assert(bits == 0);
      assert(data1.getType() >= DN_DATA);
      assert(data2.getType() >= DN_DATA);
      setType(IN_DATA);
      assert(data1.getTopology() == data2.getTopology());
      assert(data1.size() == data2.size());
      setTopology(data1.getTopology());
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = data1.size();
      size_ptr = data1.size_ptr;
      assert(data_ptr == NULL);
      data_ptr = new int[*size_ptr];
      stride = 1;
      mem_flag = true;
      access_flag = true;
      setDdeStuff(data1,data2);
    }

    void new_in(const int _topology,int &_n) {
      assert(bits == 0);
      setType(IN_DATA);
      setTopology(_topology);
      assert(size_ptr == NULL);
      //size_ptr = new int;
      //*size_ptr = _n;
      size_ptr = &_n;
      assert(data_ptr == NULL);
      data_ptr = new int[*size_ptr];
      stride = 1;
      mem_flag = true;
      access_flag = true;
      xora_ptr = NULL;
      dde_ptr = NULL;
    }

    void alloc_in() {
      //int* tmp = size_ptr;
      //size_ptr = new int;
      //*size_ptr = *tmp;
      assert(size_ptr != NULL);
      assert(data_ptr == NULL);
      data_ptr = new int[*size_ptr];
    }

    bool hasDdeStuff() const {
      if (xora_ptr == NULL) {
	assert(dde_ptr == NULL);
	return false;
      }
      else {
	assert(dde_ptr != NULL);
	return true;
      }
    }

    void setDdeStuff(const CtiData& data1,const CtiData& data2) {
      // for xora_ptr, it may exist in one or both datas -- they should be the same...
      if (data1.xora_ptr) {
	assert(data1.dde_ptr);
	// data2 should either be the same or
	assert( ((data2.xora_ptr == NULL)&&(data2.dde_ptr == NULL)) || ((data1.xora_ptr == data2.xora_ptr)&&(data1.dde_ptr == data2.dde_ptr)) );
	xora_ptr = data1.xora_ptr;
	dde_ptr = data1.dde_ptr;
      }
      else if (data2.xora_ptr) {
	assert(data2.dde_ptr);
	xora_ptr = data2.xora_ptr;
	dde_ptr = data2.dde_ptr;
      }
      else {
	assert(data1.xora_ptr == NULL);
	assert(data1.dde_ptr == NULL);
	assert(data2.xora_ptr == NULL);
	assert(data2.dde_ptr == NULL);
	assert(xora_ptr == NULL);
	assert(dde_ptr == NULL);
      }
    }

    inline string &s() {
      assert(getType() == S_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *(string*)data_ptr;
      //else
      //  return **(int**)data_ptr;
    }

    inline string s() const {
      assert(getType() == S_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *(string*)data_ptr;
      //else
      //  return **(int**)data_ptr;
    }

    inline int &i() {
      assert(getType() == I_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *(int*)data_ptr;
      //else
      //  return **(int**)data_ptr;
    }

    inline int i() const {
      assert(getType() == I_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *(int*)data_ptr;
      //else
      //  return **(int**)data_ptr;
    }

    inline double &d() {
      assert(getType() == D_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *(double*)data_ptr;
      //else
      //  return **(double**)data_ptr;
    }

    inline double d() const {
      assert(getType() == D_DATA);
      assert(data_ptr);
      //if (access_flag)
      return *(double*)data_ptr;
      //else
      //assert(0); // CI -- needed access to solver time - above return seems to work - FH
      //  return **(double**)data_ptr;
    }

    inline double &d3(const int i) {
      assert(getType() == D3_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *((double*)data_ptr+i);
      //else
      //  return *(*(double**)data_ptr+i);
    }

    inline double d3(const int i) const {
      assert(getType() == D3_DATA);
      assert(data_ptr);
      //if (access_flag)
        return *((double*)data_ptr+i);
      //else
      //  return *(*(double**)data_ptr+i);
    }

    inline double &dn(const int i) {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert(getType() == DN_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2);
      }
      else {
        assert(offset%2 == 0);
        return *(*(double**)data_ptr+i*stride/2+offset/2);
      }
    }

    inline double dn(const int i) const {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert(getType() == DN_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2);
      }
      else {
        assert(offset%2 == 0);
        return *(*(double**)data_ptr+i*stride/2+offset/2);
      }
    }

    inline double &dn3(const int i,const int j) {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert((j >= 0)&&(j < 3));
      assert(getType() == DN3_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2+j);
      }
      else {
        assert(offset%2 == 0);
        assert(*(double**)data_ptr != NULL); // TODO
        return *(*(double**)data_ptr+i*stride/2+offset/2+j);
      }
    }

    inline double dn3(const int i,const int j) const {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert((j >= 0)&&(j < 3));
      assert(getType() == DN3_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2+j);
      }
      else {
        assert(offset%2 == 0);
        assert(*(double**)data_ptr != NULL); // TODO
        return *(*(double**)data_ptr+i*stride/2+offset/2+j);
      }
    }

    inline double &dn33(const int i,const int j,const int k) {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert((j >= 0)&&(j < 3));
      assert((k >= 0)&&(k < 3));
      assert(getType() == DN33_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2+3*j+k);
      }
      else {
        assert(offset%2 == 0);
        return *(*(double**)data_ptr+i*stride/2+offset/2+3*j+k);
      }
    }

    inline double dn33(const int i,const int j,const int k) const {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert((j >= 0)&&(j < 3));
      assert((k >= 0)&&(k < 3));
      assert(getType() == DN33_DATA);
      assert(stride > 0);
      assert(stride%2 == 0);
      assert(data_ptr);
      if (access_flag) {
        return *((double*)data_ptr+i*stride/2+3*j+k);
      }
      else {
        assert(offset%2 == 0);
        return *(*(double**)data_ptr+i*stride/2+offset/2+3*j+k);
      }
    }

    inline int &in(const int i) {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert(getType() == IN_DATA);
      assert(stride > 0);
      assert(data_ptr);
      assert(size_ptr);
      if (access_flag)
        return *((int*)data_ptr+i*stride);
      else
        return *(*(int**)data_ptr+i*stride+offset);
    }

    inline int in(const int i) const {
      assert(size_ptr);
      assert((i >= 0)&&(i < *size_ptr));
      assert(getType() == IN_DATA);
      assert(stride > 0);
      assert(data_ptr);
      if (access_flag)
        return *((int*)data_ptr+i*stride);
      else
        return *(*(int**)data_ptr+i*stride+offset);
    }

    // ==================================
    // set/get functions...
    // ==================================

    void setMemFlag(const bool _mem_flag) { mem_flag = _mem_flag; }
    bool getMemFlag() const { return mem_flag; }

    void setAccessFlag(const bool _access_flag) { access_flag = _access_flag; }
    bool getAccessFlag() const { return access_flag; }

    void setFlag(const int _flag) { flag = _flag; }
    int getFlag() const { return flag;}

    /*
    bool isChtData() {
      return ( (getIndex() >= 1)&&(getUnindexedTopology() == NO_DATA) );
    }
    */

    // for field data...
    template<class T>
    void setDataPtr(T *&ptr) {
      data_ptr = (void*)&ptr;
    }
    template<class T>
    void setDataPtr(T &ptr) {
      data_ptr = (void*)&ptr;
    }

    void * getDataPtr() {
      return data_ptr;
    }

    int * getINptr() {
      // for now, it's ok to return a int pointer so long as the memory is contiguous...
      if ((getType() == IN_DATA)&&(stride == 2)) {
	assert(data_ptr);
        if ( access_flag )
          return (int*)(data_ptr);
        else
          return *(int**)data_ptr;
      }
      else {
        assert(0);
        return NULL;
      }
    }

    double * getDNptr() const {
      // for now, it's ok to return a double pointer so long as the memory is contiguous...
      if ((getType() == DN_DATA)&&(stride == 2)) {
	assert(data_ptr);
        if ( access_flag )
          return (double*)(data_ptr);
        else
          return *(double**)data_ptr;
      }
      else {
        assert(0);
        return NULL;
      }
    }

    double (* getDN3ptr() )[3] {
      if ((getType() == DN3_DATA)&&(stride == 6)) {
	assert(data_ptr);
        if ( access_flag )
          return (double(*)[3]) data_ptr;
        else
          return *(double (**)[3])data_ptr;
      }
      else {
	assert(0);
        return NULL;
      }
    }

    double (* getDN33ptr() )[3][3] {
      if ((getType() == DN33_DATA)&&(stride == 18)) {
	assert(data_ptr);
        if ( access_flag )
          return (double(*)[3][3]) data_ptr;
        else
          return *(double (**)[3][3])data_ptr;
      }
      else {
	assert(0);
        return NULL;
      }
    }

    void setStride(int _stride) {
      assert(stride == 0);
      stride = _stride;
    }

    void setOffset(int _offset) {
      assert(offset == 0);
      offset = _offset;
    }

    void setBits(const uint _bits) { bits |= _bits; }
    void clearBits(const uint _bits) { bits &= ~(_bits);}
    uint getBits() const { return bits; }
    bool checkBit(const uint _bit) { return bits & _bit; }

    void setSizePtr(int &_n) {
      assert(size_ptr == NULL);
      size_ptr = &_n;
    }
    int size() const {
      assert(size_ptr);
      return *size_ptr;
    };
    int * sizePtr() const {
      return size_ptr;
    };

    void setType(const int type);
    int getType() const;
    string getTypeAsString() const;

    void setTopology(const int topo);
    void setTopology(const int topo,const int index);
    int getTopology() const;
    string getTopologyAsString() const;
    int getUnindexedTopology() const;
    int getIndex() const;

    void setDdeStuff(int8 * &_xora,DistributedDataExchanger * &_dde);
    void initData();

    void writeData(const string& name,MPI_File& fh,MPI_Offset& offset);
    void readData(MPI_File& fh,const Header& header,const MPI_Offset& offset,const bool byte_swap = false);

    void redistReorderData(int* data);
    void redistReorderData(int8* data);
    void redistReorderData(double* data);
    void redistReorderData(double (*data)[3]);

    // these are the accessors to NData exposed...
    bool getX(double (*&x)[3]) const;
    bool getL(double *&delta) const;
    bool getGlobalIndex(int8 *&index) const;
    bool hasTets() const;
    bool getTets(int (*&noote)[4],int &nte) const;
    string getName() const;

  };

  class CtiDataProducer {
  public:
    virtual CtiDataError funcEvalCtiData(CtiData& v,const string& name,list<CtiData>& args,
                                         const bool b_eval_func) = 0;

    virtual CtiData* varEvalCtiData(const string& name) {return NULL;}

    virtual bool getSpecifierByIndex(string& name, const int index) {return false;}

    // for now this only holds the <fxn> (or <zone>:<fxn>) for argument-free fxns that return scalars
    // paranthesis aren't included
    list<string> funcEvalList;
  };

  void addCtiDataProducer(CtiDataProducer * cdp);
  void removeLastCtiDataProducer();

  struct key_lcase_equal {
    string lcs;
    key_lcase_equal(const string& s) : lcs(MiscUtils::toLowerCase(s)) {}
    bool operator()(const map<const string,CtiData>::value_type& p) const {
      return MiscUtils::toLowerCase(p.first) == lcs;
    }
  };
  inline map<const string,CtiData>::iterator find_ignore_case(map<const string,CtiData>& m,const string& s) {
    return find_if(m.begin(), m.end(), key_lcase_equal(s));
  }

  extern map<const string,CtiData> registeredDataMap; // for template class registration
  extern map<const string,CtiData> currentDataMap; // XXX stb -- is there a reason this wasnt visible before?

  void dumpCurrentData();
  void dumpRegisteredData();

  inline std::ostream& operator<< (std::ostream& os, const CtiData& data) {

    // this should only be called on rank 0...
    //assert(mpi_rank == 0);

    const int dt = data.getType();
    if (dt == I_DATA) {
      assert(data.getTopology() == 0);
      os << "I_DATA: " << data.i();
      os << " mem_flag: " << data.getMemFlag() << ", ";
    }
    else if (dt == IN_DATA) {
      os << "IN_DATA|" << data.getTopologyAsString() << ": ";
      /*
      for (int i = 0; i < min(3,data.size()); ++i) {
	cout << data.in(i) << ", ";
      }
      */
      os << "...mem_flag: " << data.getMemFlag() << " ";
    }
    else if (dt == D_DATA) {
      assert(data.getTopology() == 0);
      os << "D_DATA value: " << data.d();
      os << " mem_flag: " << data.getMemFlag() << " ";
    }
    else if (dt == D3_DATA) {
      assert(data.getTopology() == 0);
      os << "D3_DATA value: " << "[" << data.d3(0) << "," << data.d3(1) << "," << data.d3(2) << "], ";
      os << " mem_flag: " << data.getMemFlag() << " ";
    }
    else if (dt == DN_DATA) {
      os << "DN_DATA|" << data.getTopologyAsString() << ": ";
      /*
      for (int i = 0; i < min(3,data.size()); ++i) {
	cout << data.dn(i) << ", ";
      }
      */
      os << "...mem_flag: " << data.getMemFlag() << " ";
    }
    else if (dt == DN3_DATA) {
      os << "DN3_DATA|" << data.getTopologyAsString() << ": ";
      /*
      for (int i = 0; i < min(3,data.size()); ++i) {
	cout << "[" << data.dn3(i,0) << "," << data.dn3(i,1) << "," << data.dn3(i,2) << "], ";
      }
      */
      os << "...mem_flag: " << data.getMemFlag() << " ";
    }
    else if (dt == DN33_DATA) {
      os << "DN33_DATA|" << data.getTopologyAsString() << ": ";
      /*
      for (int i = 0; i < min(3,data.size()); ++i) {
	cout << "[" << data.dn33(i,0,0) << "," << data.dn33(i,0,1) << "," << data.dn33(i,0,2) << ","
                    << data.dn33(i,1,0) << "," << data.dn33(i,1,1) << "," << data.dn33(i,1,2) << ","
                    << data.dn33(i,2,0) << "," << data.dn33(i,2,1) << "," << data.dn33(i,2,2) << " ], ";
      }
      */
      os << "...mem_flag: " << data.getMemFlag() << " ";
    }
    else {
      //if (data.hash & CTI_VAR_STRING_BIT){
      //os << "string: \"" << data.str << "\"";
      assert(dt == 0);
      os << "UNDEFINED, bits: " << data.getBits() << " mem_flag: " << data.getMemFlag() << " ";
    }
    return os;

  }

  // expression parsing infrastructure...
  // getRegisteredCtiData does not carry a eval_func boolean option since the
  // call will return a ptr to the data if it exists (there is no evaluation taking place)

  CtiData * getCtiData(const string& expression, const bool eval_func = true);
  CtiData * getRegisteredCtiData(const string& expression,const bool verbose = true);
  CtiData * getUnregisteredCtiData(const string& expression, const bool eval_func = true );

  bool isRegisteredData(CtiData * data);

  // defines are supported in a key/value tring map...
  //void addDefine(const string& key,const string& value);

  // ==================================
  // interacting ith registered data...
  // ==================================

  void setRegisteredCvDnAndDn3Names(vector<string>& nameVec);
  void setRegisteredDAndINames(vector<string>& nameVec);
  void setRegisteredCvAndNoDNNames(vector<string>& nameVec);
  void setRegisteredCvAndNoDN3Names(vector<string>& nameVec);
  void setRegisteredBfDNNames(vector<string>& nameVec);
  void setRegisteredBfDN3Names(vector<string>& nameVec);
  void setRegisteredCvBfAndNoDNNames(vector<string>& nameVec);
  void setRegisteredCvBfAndNoDN3Names(vector<string>& nameVec);
  void setRegisteredFuncEvalNames(vector<string>& nameVec);

  // ===========================================
  // int value I
  // ===========================================

  void _registerData(int& val,const string& name,const uint rw_bits);

  // ==================================
  // double value D
  // ==================================

  void _registerData(double& val,const string& name,const uint rw_bits);

  // ==================================
  // double[3] value D
  // ==================================

  void _registerData(double val[3],const string& name,const uint rw_bits);

  // ===================================================================================
  // double scalar DN registration
  // ===================================================================================

  void _registerData(double *&val,const string& name,const int topology,const uint rw_bits,int &n);
  void _registerData(double *&val,const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde);

  // trying to eliminate the need for topology for these registrations...
  // cannot do this because of conflict with D3_DATA registration - not completely understood...
  //void _registerData(double *&val,const string& name,const uint rw_bits,int &n);

  // ===================================================================================
  // double vector DN3 registration
  // ===================================================================================

  void _registerData(double (*&val)[3],const string& name,const int topology,const uint rw_bits,int &n);
  void _registerData(double (*&val)[3],const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde);

  // trying to eliminate the need for topology for these registrations...
  //void _registerData(double (*&val)[3],const string& name,const uint rw_bits,int &n);

  // ===================================================================================
  // int scalar IN registration
  // ===================================================================================

  void _registerData(int *&val,const string& name,const int topology,const uint rw_bits,int &n);
  void _registerData(int *&val,const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde);

  // ===================================================================================
  // versions of registration where memory is managed by the registration...
  // allows registration of variables from the input file, like prepro's REGISTER_CV_DN...
  // ===================================================================================

  void _registerData(const string& name,const int datatype,const uint rw_bits);
  void _registerData(const string& name,const int topology,const int datatype,const uint rw_bits,int &n);
  void _registerData(const string& name,const int topology,const int datatype,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde);
  void _initData();

  // ===================================================================================
  // associates an index or name to an "n" datatype...
  // ===================================================================================

  void _registerGlobalIndex(int8 * &index,int &n);
  void _registerName(const string& name,int &n);
  void _registerNoote(int (*&noote)[4],int &nte,int &n);

  // ===================================================================================
  // versions of registration where the memory is part of struct/class...
  // ===================================================================================

  template<class T>
  void _registerData(T *&state,int &val,const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array int \"" << name << "\", struct_size " << sizeof(T) << " without io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset(&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(IN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

  }

  template<class T>
  void _registerData(T *&state,int &val,const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array int \"" << name << "\", struct_size " << sizeof(T) << " with io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset(&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(IN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

  }

  template<class T>
  void _registerData(T *&state,double &val,const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array double \"" << name << "\", struct_size " << sizeof(T) << " without io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset((int*)&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(DN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

  }

  template<class T>
  void _registerData(T *&state,double &val,const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array double \"" << name << "\", struct_size " << sizeof(T) << " with io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset((int*)&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(DN_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

  }

  template<class T>
  void _registerData(T *&state,double (&val)[3],const string& name,const int topology,const uint rw_bits,int &n) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array double[3] \"" << name << "\", struct_size " << sizeof(T) << " without io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset((int*)&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(DN3_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // ensure that dde is provided whenever i/o is needed...

    assert(!( (rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA) ));

  }

  template<class T>
  void _registerData(T *&state,double (&val)[3],const string& name,const int topology,const uint rw_bits,int &n,int8 * &xora,DistributedDataExchanger * &dde) {

    if (mpi_rank == 0) cout << " > CtiRegister: struct array double[3] \"" << name << "\", struct_size " << sizeof(T) << " with io" << endl;

    // stride check...

    if (sizeof(T)%sizeof(double) != 0) {
      CERR("stride alignment problem: " << sizeof(T) << " " << sizeof(double));
    }

    // map it...

    pair<const string,CtiData> key_value_pair = pair<const string,CtiData>(name,CtiData());
    pair<map<const string,CtiData>::iterator,bool> return_pair = registeredDataMap.insert(key_value_pair);

    // uniquesness check...

    if (!return_pair.second) {
      CERR("duplicate name problem \"" << name << "\"");
    }
    assert(key_value_pair.second.empty());

    // build the crd...

    CtiData * data = &(return_pair.first->second);
    data->setDataPtr(state);
    data->setStride(sizeof(T)/sizeof(int));
    data->setOffset((int*)&val-(int*)state);
    data->setBits(rw_bits);
    data->setType(DN3_DATA);
    data->setTopology(topology);
    data->setSizePtr(n);

    // dde consistency check...

    assert((rw_bits & READ_DATA) || (rw_bits & WRITE_DATA) || (rw_bits & CAN_WRITE_DATA));
    //assert(xora);
    //assert(dde);
    //assert(n == dde->my_nda);
    data->setDdeStuff(xora,dde);

  }

  string getFirstSpecifier(const string var,const bool filter_protected=true);  // get stuff before ":"
  string getSecondSpecifier(const string var);  // get stuff after ":"

  // data flag

  int getDataFlag(const string& name);
  void setDataFlag(const string& name,const int val);
  bool checkDataFlag(const string& name);
  void clearAllDataFlags();

  // RESULT io

  void readData(const string& filename);
  void readLpData(const string& filename);
  void writeData(const string& filename);
  void writeData(MPI_File &fh, MPI_Offset &offset);

  // STATS...

  class Stats {
  public:
    string name;
    int datatype;
    // pointers into the registered data map...
    CtiData * source;
    CtiData * wgt;
    CtiData * avg;
    CtiData * rms;
    CtiData * rey;
    Stats(const string& _name,const int _datatype) : name(_name), datatype(_datatype) {
      source = NULL;
      wgt = NULL;
      avg = NULL;
      rms = NULL;
      rey = NULL;
    }
  };

  extern vector<Stats> statsVec;
  void registerStats(const string& vname,CtiData * data,const bool b_init);
  void _initStats();
  void resetStats();
  void updateStats(const double wgt,const bool verbose=false);
  void clearStatsBits(const uint _bits);

  // time step actions...

  void clearCurrentData();

}

#endif
