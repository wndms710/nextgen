#ifndef FLUENT_PARSER_HPP
#define FLUENT_PARSER_HPP

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "ByteSwap.hpp"
#include "MiscUtils.hpp"

using namespace std;

/*
This is a base class that can be used for parsing Fluent formatted files. It includes ascii and binary chunked parsing, and is designed for both serial or parallel use (controlled via compilation flags). Various file (i.e., .cas, .dat, .msh) readers can be derived from this.

The main difference between this and a generic chunked file reader is that it comprehends "levels", i.e., depth of nested parentheses, which the Fluent format uses to delineate blocks of data.
*/

#define FLUENT_BUFFER_SIZE     65536
#define FLUENT_TOKEN_SIZE      65536

class FluentParser {

private:

  int pos, max_pos;
  char buf[FLUENT_BUFFER_SIZE+FLUENT_TOKEN_SIZE];

protected:

  FILE * fp;
  int level;
  bool byteSwap;

public:

  FluentParser() {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    level = 0;
    fp = NULL;
    byteSwap = false;
  }

  FluentParser(const string& filename) {
    init(filename);
  }

  ~FluentParser() {
    if (fp != NULL) fclose(fp);
  }

  void init(const string& filename) {
    fp = NULL;  // properly intialize to NULL for non-read ranks
    COUT2("FluentParser::init(" << filename << ")");
    // read from rank0 only and bcast buffer for parallel I/O
    if (mpi_rank == 0) {
      const int file_err = MiscUtils::openFile(&fp,filename,"rb");
      if (file_err != 0) throw(file_err);
    }

    // reset tokenizing stuff...
    pos = max_pos = 0;
    level = 0;
    byteSwap = false;

  }

  void setByteSwap(const bool _byteSwap) {
    byteSwap = _byteSwap;
  }

  inline int getLevel() const {return level;}

  void rewindFile() {
    if (mpi_rank == 0) rewind(fp);

    // also need to rewind position stuff; should trigger re-buffering from start
    pos = max_pos = level = 0;
  }

  /*
  ASCII parsing of file
  */

  // only rank0 reads the file into chunked buffer, however everybody does the same parsing once the buffer has been broadcast
  int getNextToken(char *& token) {

    token = NULL;
    // int token_pos = 0;
    int token_level = level;
    int quote_mode = 0;

    while (1) {

      if (pos >= max_pos) {

        if (token) {
          // need to shift buffer but token isn't done yet, so shift what is in token to start of buf and reset token location
          pos = max_pos - int(token-buf);  // bytes of token already read
          if (token != buf) {
            memmove(buf,token,pos);
            token = buf; // reset to start...
          }
        }
        else {
          // just reset pos to start of buffer
          pos = 0;
        }

        if (mpi_rank == 0) max_pos = pos + fread(buf+pos, sizeof(char), FLUENT_BUFFER_SIZE-pos, fp);

        MPI_Bcast(buf,FLUENT_BUFFER_SIZE,MPI_CHAR,0,mpi_comm);
        MPI_Bcast(&max_pos,1,MPI_INT,0,mpi_comm);

        if (max_pos == pos) {
          buf[pos] = '\0';
          cerr << "Error at end of file: \"" << buf << "\"" << endl;
          assert(0);
        }
      }

      const char c = buf[pos++];

      if (c == '"') {
        if (quote_mode == 1) break;

        token = buf+pos; // first char after quote
        quote_mode = 1;
      }
      else if (quote_mode == 0) {
        if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) || (c == '\r') ) {
          // c is a whitespace character - this represents either the termination
          // of a token, or some space at the start of the next token...
          if (token) break;
        }
        else if (c == '(') {
          ++level;
          if (token) break;
          ++token_level;
        }
        else if (c == ')') {
          --level;
          if (token) break;
          --token_level;
        }
        else if (!token) {
          // any other character is the start of the token...
          token = buf+pos-1;
        }
      }
    }

    // terminate string and return...
    buf[pos-1] = '\0';

    return (token_level);
  }

  // move through ASCII tokens until specific end of level criteria are hit
  bool advanceToLevel(const int target_level) {
    int quote_mode = 0;

    while (1) {

      if (pos >= max_pos) {

        if (mpi_rank == 0) max_pos = fread(buf, sizeof(char), FLUENT_BUFFER_SIZE, fp);

        MPI_Bcast(buf,FLUENT_BUFFER_SIZE,MPI_CHAR,0,mpi_comm);
        MPI_Bcast(&max_pos,1,MPI_INT,0,mpi_comm);

        if (max_pos == 0) return false;
        pos = 0;
      }

      char c = buf[pos++];

      if (c == '"') {
        quote_mode = 1 - quote_mode;
      }
      else if (quote_mode == 0) {
        if (c == '(') {
          ++level;
          if (level == target_level) return true;
        }
        else if (c == ')') {
          --level;
        }
      }
    }
  }

  /*
  Binary parsing methods
  */

  template<typename T_To, typename T_From>
  T_To getNextBinary() {
    T_From d;
    // we need to read sizeof(T_From) bytes from the current buf+pos. If the remaining
    // bytes are not sufficient to read, then realign before reading...

    if (pos+int(sizeof(T_From)) > max_pos) {
      assert(pos > 0);

      // only rank0 computes shifted buffer and reads; then send
      if ( mpi_rank == 0 ) {
        const int n = max_pos-pos;
        if (n > 0) memmove(buf,buf+(pos),n);
        max_pos = n + fread(buf+n, sizeof(char), FLUENT_BUFFER_SIZE-n, fp);
      }
      MPI_Bcast(buf,FLUENT_BUFFER_SIZE,MPI_CHAR,0,mpi_comm);
      MPI_Bcast(&max_pos,1,MPI_INT,0,mpi_comm);

      pos = 0;
    }
    assert(pos+int(sizeof(T_From)) <= max_pos);
    memcpy(&d,buf+pos,sizeof(T_From));
    if (byteSwap) d = ByteSwap::byteSwap(d);
    pos += sizeof(T_From);

    return T_To(d);  // cast to the appropriate format
  };

  void readBinaryInt(int& d) {
    d = getNextBinary<int,int>();
  }

  /*
  methods to convert char * to various types
  */

  inline int asInt(const char * token) {
    return atoi(token);
  }

  inline double asDouble(const char * token) {
    double d;
    sscanf(token, "%lf", &d);
    return (d);
  }

  inline string asString(const char * token) {
    return string(token);
  }

  inline int asHex(const char * token) {
    int i;
    sscanf(token, "%x", &i);
    return (i);
  }

  inline int getNextTokenAsInt(const int level_check) {
    char * tmp_token = NULL;
    const int level = getNextToken(tmp_token); assert(level == level_check);
    return asInt(tmp_token);
  }

  inline double getNextTokenAsDouble(const int level_check) {
    char * tmp_token = NULL;
    const int level = getNextToken(tmp_token); assert(level == level_check);
    return asDouble(tmp_token);
  }

  inline string getNextTokenAsString(const int level_check) {
    char * tmp_token = NULL;
    const int level = getNextToken(tmp_token); assert(level == level_check);
    return asString(tmp_token);
  }

  inline int getNextTokenAsHex(const int level_check) {
    char * tmp_token = NULL;
    const int level = getNextToken(tmp_token); assert(level == level_check);
    return asHex(tmp_token);
  }

};

#endif
