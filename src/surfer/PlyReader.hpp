#ifndef PLY_READER_HPP
#define PLY_READER_HPP

#include <math.h>
#include <assert.h>
#include "ByteSwap.hpp"

using namespace std;

#define PLY_BUFFER_SIZE     65536
#define PLY_TOKEN_SIZE      65536

class PlyProperty {
private:
  int * pos;
  int * max_pos;
  char * buf;
  FILE * fp;
  bool b_ascii;
  bool byteSwap;

  bool b_isList;
  bool b_isInt;
  string var_type;
  string list_var_type;
  string name;
  int (PlyProperty::*readListInt)();
  int (PlyProperty::*readInt)();
  double (PlyProperty::*readDouble)();
public:
  PlyProperty(const bool _b_ascii,int& _pos,int& _max_pos,char (&_buf)[PLY_BUFFER_SIZE+PLY_TOKEN_SIZE],FILE * _fp,const bool _byteSwap) {
    pos = &_pos;
    max_pos = &_max_pos;
    b_ascii = _b_ascii;
    byteSwap = _byteSwap;
    buf = _buf;
    fp = _fp;

    b_isList = false;
    b_isInt = false;
    var_type = "";
    list_var_type = "";
    name = "";
    readInt = NULL;
    readListInt = NULL;
    readDouble = NULL;
  };

  ~PlyProperty() {
    readInt = NULL;
    readListInt = NULL;
    readDouble = NULL;
    fp = NULL;
  };

  void define(const string _list_var_type,const string _var_type,const string _name) {
    b_isList = true;
    var_type = _var_type;
    list_var_type = _list_var_type;
    name = _name;

    if (!b_ascii) {
      if (list_var_type == "char" || list_var_type == "int8") readListInt = &PlyProperty::readBinary<int,char>;
      else if (list_var_type == "uchar" || list_var_type == "uint8") readListInt = &PlyProperty::readBinary<int,unsigned char>;
      else if (list_var_type == "short" || list_var_type == "int16") readListInt = &PlyProperty::readBinary<int,short int>;
      else if (list_var_type == "ushort" || list_var_type == "uint16") readListInt = &PlyProperty::readBinary<int,unsigned short int>;
      else if (list_var_type == "int" || list_var_type == "int32") readListInt = &PlyProperty::readBinary<int,int>;
      else if (list_var_type == "uint" || list_var_type == "uint32") readListInt = &PlyProperty::readBinary<int,unsigned int>;

      if (var_type == "char" || var_type == "int8") {
        readInt = &PlyProperty::readBinary<int,char>;
        b_isInt = true;
      }
      else if (var_type == "uchar" || var_type == "uint8") {
        readInt = &PlyProperty::readBinary<int,unsigned char>;
        b_isInt = true;
      }
      else if (var_type == "short" || var_type == "int16") {
        readInt = &PlyProperty::readBinary<int,short int>;
        b_isInt = true;
      }
      else if (var_type == "ushort" || var_type == "uint16") {
        readInt = &PlyProperty::readBinary<int,unsigned short int>;
        b_isInt = true;
      }
      else if (var_type == "int" || var_type == "int32") {
        readInt = &PlyProperty::readBinary<int,int>;
        b_isInt = true;
      }
      else if (var_type == "uint" || var_type == "uint32") {
        readInt = &PlyProperty::readBinary<int,unsigned int>;
        b_isInt = true;
      }
      else if (var_type == "float" || var_type == "float32") readDouble = &PlyProperty::readBinary<double,float>;
      else if (var_type == "double" || var_type == "float64") readDouble = &PlyProperty::readBinary<double,double>;
    }
    else {
      readListInt = &PlyProperty::getNextTokenAsInt;

      if (var_type == "char" || var_type == "int8" || var_type == "uchar" || var_type == "uint8" || var_type == "short" || var_type == "int16" || var_type == "ushort" || var_type == "uint16" || var_type == "int" || var_type == "int32" || var_type == "uint" || var_type == "uint32") {
          readInt = &PlyProperty::getNextTokenAsInt;
          b_isInt = true;
      }
      else if (var_type == "float" || var_type == "float32" || var_type == "double" || var_type == "float64") {
        readDouble = &PlyProperty::getNextTokenAsDouble;
      }
    }
  };
  void define(const string _var_type,const string _name) {
    var_type = _var_type;
    name = _name;

    if (!b_ascii) {
      if (var_type == "char" || var_type == "int8") {
        readInt = &PlyProperty::readBinary<int,char>;
        b_isInt = true;
      }
      else if (var_type == "uchar" || var_type == "uint8") {
        readInt = &PlyProperty::readBinary<int,unsigned char>;
        b_isInt = true;
      }
      else if (var_type == "short" || var_type == "int16") {
        readInt = &PlyProperty::readBinary<int,short int>;
        b_isInt = true;
      }
      else if (var_type == "ushort" || var_type == "uint16") {
        readInt = &PlyProperty::readBinary<int,unsigned short int>;
        b_isInt = true;
      }
      else if (var_type == "int" || var_type == "int32") {
        readInt = &PlyProperty::readBinary<int,int>;
        b_isInt = true;
      }
      else if (var_type == "uint" || var_type == "uint32") {
        readInt = &PlyProperty::readBinary<int,unsigned int>;
        b_isInt = true;
      }
      else if (var_type == "float" || var_type == "float32") readDouble = &PlyProperty::readBinary<double,float>;
      else if (var_type == "double" || var_type == "float64") readDouble = &PlyProperty::readBinary<double,double>;
    }
    else {
      if (var_type == "char" || var_type == "int8" || var_type == "uchar" || var_type == "uint8" || var_type == "short" || var_type == "int16" || var_type == "ushort" || var_type == "uint16" || var_type == "int" || var_type == "int32" || var_type == "uint" || var_type == "uint32") {
          readInt = &PlyProperty::getNextTokenAsInt;
          b_isInt = true;
      }
      else if (var_type == "float" || var_type == "float32" || var_type == "double" || var_type == "float64") {
        readDouble = &PlyProperty::getNextTokenAsDouble;
      }
    }
  };

  inline int getListInt() {return (*this.*readListInt)();};
  inline int getInt() {return (*this.*readInt)();};
  inline double getDouble() {return (*this.*readDouble)();};
  bool isInt() const {return b_isInt;};
  bool isList() const {return b_isList;};
  string getName() const {return name;};

  template<typename T_To, typename T_From>
  T_To readBinary() {
    T_From d;
    // we need to read sizeof(T_From) bytes from the current buf+pos. If the remaining
    // bytes are not sufficient to read, then realign before reading...

    if (*pos+int(sizeof(T_From)) > *max_pos) {
      assert(*pos > 0);
      const int n = *max_pos-*pos;
      if (n > 0)
        memmove(buf,buf+(*pos),n);

      *max_pos = n + fread(buf+n, sizeof(char), PLY_BUFFER_SIZE, fp);
      *pos = 0;
    }
    assert(*pos+int(sizeof(T_From)) <= *max_pos);
    memcpy(&d,buf+(*pos),sizeof(T_From));
    if (byteSwap) d = ByteSwap::byteSwap(d);
    *pos += sizeof(T_From);

    return T_To(d);
  };

  char * getNextPlyToken() {

    char * token = NULL;
    while (1) {

      if (*pos >= *max_pos) {

        assert(*pos == *max_pos);
        if (token) {
          // if the token is active, but not completed, then we need to shift the token part of
          // the current buf (i.e. the end) to the start of the buf, and read the next part in...
          //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
          *pos = *max_pos - int(token-buf);
          if (token != buf) {
            memmove(buf,token,*pos);
            token = buf; // reset to start...
          }
        }
        else {
          *pos = 0;
        }

        *max_pos = *pos + fread(buf+(*pos), sizeof(char), PLY_BUFFER_SIZE, fp);
        if (*max_pos == *pos) {
          buf[*pos] = '\0';
          return token;
        }
      }

      const char c = buf[(*pos)++];

      if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) ) {
        // c is a whitespace character - this represents either the termination
        // of a token, or some space at the start of the next token...
        if (token) break;
      }
      else if (!token) {
        // any other character is the start of the token...
        token = buf+(*pos)-1;
      }

    }

    // terminate string and return...
    buf[*pos-1] = '\0';
    return token;
  }

  int getNextTokenAsInt() {
    char * token = getNextPlyToken();
    return atoi(token);
  }

  string getNextTokenAsString() {
    char * token = getNextPlyToken();
    return string(token);
  }

  double getNextTokenAsDouble() {
    char * token = getNextPlyToken();
    return atof(token);
  }

};

struct PlyElement {
private:
  int count;
  string name;
public:
  vector<PlyProperty> propertiesVec;
  PlyElement() {
    count = 0;
    name = "";
  };
  PlyElement(const string _name,const int _count) {
    name = _name;
    count = _count;
  };
  ~PlyElement() {
    propertiesVec.clear();
  };

  void setCount(const int _count) {count = _count;};
  int getCount() const {return count;};

  void setName(const int _name) {name = _name;};
  string getName() const {return name;};
};

class PlyReader {
private:

  int pos, max_pos;
  bool b_ascii;
  bool byteSwap;
  char buf[PLY_BUFFER_SIZE+PLY_TOKEN_SIZE];

protected:

  FILE * fp;

public:

  vector<PlyElement> elementsVec;

  PlyReader() {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    byteSwap = false;
    b_ascii = true;
    fp = NULL;
  }

  PlyReader(const string& filename) {
    init(filename);
  }

  ~PlyReader() {
    if (fp != NULL) fclose(fp);
    elementsVec.clear();
  }

  void init(const string& filename) {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    byteSwap = false;
    b_ascii = true;

    const int file_err = MiscUtils::openFile(&fp,filename,"rb");
    if (file_err != 0) throw(file_err);

    // fp = fopen(filename.c_str(),"rb");
    // if (fp == NULL) {
    //   cerr << "Error: cannot open file: "<< filename << endl;
    //   throw(20);
    // }
  }

  void rewindFile() {
    rewind(fp);
  }

  void setBinary() {
    b_ascii = false;
  }

  void setByteSwap() {
    byteSwap = true;
  }

  void addElement(const string name,const int count) {
    elementsVec.push_back(PlyElement(name,count));
  }

  void addProperty(const string var_type,const string name) {
    elementsVec.back().propertiesVec.push_back(PlyProperty(b_ascii,pos,max_pos,buf,fp,byteSwap));
    elementsVec.back().propertiesVec.back().define(var_type,name);
  }

  void addProperty(const string list_var_type,const string var_type,const string name) {
    elementsVec.back().propertiesVec.push_back(PlyProperty(b_ascii,pos,max_pos,buf,fp,byteSwap));
    elementsVec.back().propertiesVec.back().define(list_var_type,var_type,name);
  }

  char * getNextPlyToken() {

    char * token = NULL;
    while (1) {

      if (pos >= max_pos) {

        assert(pos == max_pos);
        if (token) {
          // if the token is active, but not completed, then we need to shift the token part of
          // the current buf (i.e. the end) to the start of the buf, and read the next part in...
          //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
          pos = max_pos - int(token-buf);
          if (token != buf) {
            memmove(buf,token,pos);
            token = buf; // reset to start...
          }
        }
        else {
          pos = 0;
        }

        max_pos = pos + fread(buf+pos, sizeof(char), PLY_BUFFER_SIZE, fp);
        if (max_pos == pos) {
          buf[pos] = '\0';
          return token;
        }
      }

      const char c = buf[pos++];

      if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) ) {
        // c is a whitespace character - this represents either the termination
        // of a token, or some space at the start of the next token...
        if (token) break;
      }
      else if (!token) {
        // any other character is the start of the token...
        token = buf+pos-1;
      }

    }

    // terminate string and return...
    buf[pos-1] = '\0';
    return token;
  }

  void goToNextLine() {

    char * token = NULL;
    while (1) {

      if (pos >= max_pos) {

        assert(pos == max_pos);
        if (token) {
          // if the token is active, but not completed, then we need to shift the token part of
          // the current buf (i.e. the end) to the start of the buf, and read the next part in...
          //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
          pos = max_pos - int(token-buf);
          if (token != buf) {
            memmove(buf,token,pos);
            token = buf; // reset to start...
          }
        }
        else {
          pos = 0;
        }

        max_pos = pos + fread(buf+pos, sizeof(char), PLY_BUFFER_SIZE, fp);
        if (max_pos == pos) {
          buf[pos] = '\0';
          return;
        }
      }

      const char c = buf[pos++];

      if ( (c == '\n') ) {
        // c is a whitespace character - this represents either the termination
        // of a token, or some space at the start of the next token...
        if (token) break;
      }
      else if (!token) {
        // any other character is the start of the token...
        token = buf+pos-1;
      }

    }

    // terminate string and return...
    buf[pos-1] = '\0';
    return;
  }

  int getNextTokenAsInt() {
    char * token = getNextPlyToken();
    return atoi(token);
  }

  string getNextTokenAsString() {
    char * token = getNextPlyToken();
    return string(token);
  }

  double getNextTokenAsDouble() {
    char * token = getNextPlyToken();
    return atof(token);
  }

};

#endif
