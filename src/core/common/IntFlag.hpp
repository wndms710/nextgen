#ifndef _INT_FLAG_HPP_
#define _INT_FLAG_HPP_

#include <assert.h>
#include <set>
#include <string.h>
#include <iostream>

#define FOR_IFL for (int ifl = 0; ifl < length; ++ifl)

// flag base class
class IntFlag {
private:
  int length;
  int max_length;
  int * flag;

public:

  IntFlag() {
    flag = NULL;
    length = max_length = 0;
  }

  IntFlag(int _length) {
    flag = NULL;
    length = max_length = 0;
    setLength(_length);
  }

  IntFlag(const IntFlag& other) {
    // custom copy constructor to properly re-allocate memory for copy;
    // otherwise pointing to memory that will be cleared
    flag = NULL;
    length = max_length = 0;
    setLength(other.length);
    FOR_IFL flag[ifl] = other.flag[ifl];
    // memcpy(flag,other.flag,sizeof(int)*other.length);
  }

  ~IntFlag() {
    if (flag != NULL) delete[] flag;
  }

  inline int& operator[] (int ifl) const {
    if (flag == NULL) {
      std::cout << "Error: IntFlag operator[]: flag == NULL" << std::endl;
      throw(0);
    }
    if (!((ifl >= 0) && (ifl < length))) {
      std::cout << "Error: IntFlag operator[] index: ifl: " << ifl << " out of range. length: " << length << std::endl;
      throw(0);
    }
    //if (!(ifl >= 0 && ifl < length)) CWARN("flag accessor out-of-bounds");
    //assert(flag);
    //if (!((ifl >= 0) && (ifl < length)))
    //  std::cout << "ifl: " << ifl << " length: " << length << std::endl;
    //assert((ifl >= 0) && (ifl < length));
    return flag[ifl];
  }

  bool isNull() const { return flag == NULL; }

  void setLength(int _length) {
    if (_length <= 0) return;  // invalid length specified

    if ((_length == length) && (flag != NULL) ) return;  // same length

    length = _length;
    if (length <= max_length) return;  // current allocation (max_length) is sufficient

    // flag is either null or length > max_length
    if ( (flag != NULL) && (length > max_length) )
      delete[] flag;

    max_length = length;
    flag = new int[length];
  }

  // eventually deprecate getLength() for size
  int getLength() const {
    return length;
  }

  int getMaxLength() const {
    return max_length;
  }

  int size() const {
    return length;
  }

  void resize(int _length) {

    if ( (flag == NULL) || (_length <= max_length)) {
      setLength(_length);  // if new or smaller than available memory, simply change length of index w/o reallocating
      return;
    }

    // get here means _length > max_length
    int * flag_new = new int[_length];
    memcpy(flag_new,flag,sizeof(int)*length);

    length = max_length = _length;

    delete[] flag; flag = flag_new;
  }

  void ensureSize(const int _max_length) {
    // increases the size of the flag's memory without increasing its size...
    if (_max_length > max_length) {
      assert(length <= max_length);
      max_length = _max_length;
      if (flag) {
        int * flag_new = new int[max_length];
        memcpy(flag_new,flag,sizeof(int)*length);
        delete[] flag; flag = flag_new;
      }
      else {
        assert(length == 0); // must be true
        flag = new int[max_length];
      }
    }
  }

  void clear() {
    length = max_length = 0;
    if (flag != NULL) {
      delete[] flag;
      flag = NULL;
    }
  }

  // counts non-zeros
  int count() const {
    if ( (flag == NULL) ) return 0;
    int count = 0;
    FOR_IFL {
      if (flag[ifl] != 0) count++;
    }
    return count;
  }

  int countUnique() const {
    if ( (flag == NULL) ) return 0;
    std::set<int> uniqueVals;
    FOR_IFL uniqueVals.insert(flag[ifl]);
    return uniqueVals.size();
  }

  int countPositive() const {
    if ( (flag == NULL) ) return 0;
    int count = 0;
    FOR_IFL {
      if (flag[ifl] > 0) count++;
    }
    return count;
  }

  int countNegative() const {
    if ( (flag == NULL) ) return 0;
    int count = 0;
    FOR_IFL {
      if (flag[ifl] < 0) count++;
    }
    return count;
  }

  int countZero() const {
    if ( (flag == NULL) ) return 0;
    int count = 0;
    FOR_IFL {
      if (flag[ifl] == 0) count++;
    }
    return count;
  }

  int countEqualTo(const int val) const {
    if ( (flag == NULL) ) return 0;
    int count = 0;
    FOR_IFL {
      if (flag[ifl] == val) count++;
    }
    return count;
  }

  void set(int ifl,int value) {
    if (length == 0) return;

    assert( (flag != NULL) );
    assert((ifl >= 0)&&(ifl < length));
    flag[ifl] = value;
  }

  void setAll(int value) {
    if (length == 0) return;
    assert( (flag != NULL) );
    FOR_IFL {
      flag[ifl] = value;
    }
  }

  int get(int ifl) const {
    assert( (flag != NULL) );
    assert((ifl >= 0)&&(ifl < length));
    return flag[ifl];
  }
};

#undef FOR_IFL

#endif
