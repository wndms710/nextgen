#ifndef ASCII_READER_HPP
#define ASCII_READER_HPP

#include <math.h>
#include <assert.h>
#include "ByteSwap.hpp"

using namespace std;

#define ASCII_BUFFER_SIZE     65535
#define ASCII_TOKEN_SIZE      65535

/*
  Chunked buffer parsing of ASCII files
  Allows user to:
   - request next value as specific type
   - detect end of line
   - detect end of file

   Currently end-of-line and space (token separators) characters are hard coded
   TODO: allow user to pass in these characters on construction of the reader
*/

class AsciiReader {
protected:
  FILE * fp;
  bool b_eof;        // is position at end of file
  bool b_eol;        // is position at end of line
  bool b_max_is_eof; // buffer max is end of file
  int pos, max_pos;  // current position, buffer's max position

  char buf[ASCII_BUFFER_SIZE+ASCII_TOKEN_SIZE];

public:

  AsciiReader() {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    b_eof = b_eol = false;
    b_max_is_eof = false;
    fp = NULL;
  }

  AsciiReader(const string& filename) {
    init(filename);
  }

  ~AsciiReader() {
    if (fp != NULL) fclose(fp);
  }

  bool atEof() const {
    return b_eof;
  }

  void init(const string& filename) {
    pos = max_pos = 0;
    b_eof = b_eol = false;
    b_max_is_eof = false;

    const int file_err = MiscUtils::openFile(&fp,filename,"rb");
    if (file_err != 0) throw(file_err);
  }

  void rewindFile() {
    // return position pointer to beginning of file
    rewind(fp);
    pos = max_pos = 0;
  }

  char * getNextToken(const bool b_eol_stop = false) {
    // returns next token as char array

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

        max_pos = pos + fread(buf+pos, sizeof(char), ASCII_BUFFER_SIZE, fp);
        if (feof(fp)) b_max_is_eof = true;

        if (max_pos == pos) {
          buf[pos] = '\0';
          if (b_max_is_eof) {
            b_eof = true;
          }
          return token;
        }
      }

      const char c = buf[pos++];
      // cout << "c: " << c << " " << (int)c << endl;  // super verbose debugging

      if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) ) {
        // if c is a whitespace character (valid values above)
        // this represents either the termination of a token,
        // or some space at the start of the next token...

        if ( (c=='\n') ) {
          // WARNING: including 'carriage return c==13' above causes issues with
          // proper line ending recognition; seems as though CR and NL (new line) appear
          // in succession causing issues
          b_eol = true;  // only valid end-of-line character
          if (b_eol_stop) break;
        }

        if (token) break;
      }
      else if (!token) {
        // any other character is the start of the token...
        b_eol = false;
        token = buf+pos-1;
      }

    }

    // terminate string and return...
    buf[pos-1] = '\0';

    if ((pos == max_pos) && b_max_is_eof) b_eof = true;  // if we just parsed last token in file
    return token;
  }

  void goToNextLine(const bool b_force=false) {
    // when a parameter on the line has been parsed, this acts as goToEndOfLine:
    // after parsing the final token a '\n' is found and b_eol=true and pointer is at
    // start of next line.

    // when requested again, we assume user doesn't know whether previous token parsed
    // was last in previous line. So if b_eol, we simply return - we have moved to the
    // next line but user didn't realize it.

    // the confusion is when the user calls this and they are sitting at the start of
    // the next line - how can we tell they want to really move to the NEXT line vs be
    // on the current one? Here we use the b_force to indicate skipToNextLine...but
    // usage requires user to comprehend that they are at the start of the line (not ideal)

    //TODO: think of a better way to capture line position so single goToNextLine behaves
    // intuitively regardless

    if (b_force) b_eol = false;
    if (b_eol) return;

    char * token = NULL;
    while ( !(b_eol || b_eof) ) {
      token = getNextToken(true);  // true forces parsing stop after eol found
    }
  }

  int getNextTokenAsInt() {
    char * token = getNextToken();
    return atoi(token);
  }

  string getNextTokenAsString() {
    char * token = getNextToken();
    return string(token);
  }

  double getNextTokenAsDouble() {
    char * token = getNextToken();
    return atof(token);
  }

};

#undef ASCII_BUFFER_SIZE
#undef ASCII_TOKEN_SIZE

#endif
