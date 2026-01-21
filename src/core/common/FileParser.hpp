#ifndef FILEPARSER_HPP
#define FILEPARSER_HPP

// ======================================================================
// FileParser.hpp
// 
// this class tokenizes a file one token at a time based on 
// user-defined delimiters and whitespace. Each token gets returned
// in a passed string along with the token type.  
//
// Sample usage:
// 
//  FileParser fp(filename," \t","()"); // filename, whitespace, braces
//  string token;
//  while (1) {
//    FileParserTokenType ptt = fp.getNextToken(token);
//    if (ptt == FILEPARSER_NO_TOKEN)
//      break;
//    else if (ptt == FILEPARSER_DELIMITER_TOKEN)
//      cout << "got delimiter: " << token << endl;
//   }
//
// History:
//  author: Frank Ham - March 2013
// ======================================================================

// ======================================================================
// Copyright (c) 2013 Cascade Technologies Inc.
// All Rights Reserved. www.cascadetechnologies.com
//
// This source is proprietary software and cannot be used, redistributed,
// or modified except under the terms of Cascade's Software License 
// Agreement. 
//
// THIS SOURCE CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY 
// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
// ======================================================================

#include "Common.hpp"
#include "Params.hpp"

enum FileParserTokenType {
  FILEPARSER_NO_TOKEN,
  FILEPARSER_STANDARD_TOKEN,
  FILEPARSER_DELIMITER_TOKEN
};

class FileParser {
  
private:
  
  string whitespace;
  string delimiters;
  string whitespace_and_delimiters;
  int linenum;
  string line;
  string::size_type pos;
  ifstream ifile;
  
public:
  
  FileParser(const string& filename,const string& whitespace,const string& delimiters) {
    
    ifile.open(filename.c_str());
    
    if (!ifile.is_open()) {
      cerr << "Error: FileParser could not open file \"" << filename << "\"" << endl;
      throw(-1);
    }
    
    this->whitespace = whitespace;
    this->delimiters = delimiters;
    whitespace_and_delimiters = whitespace + delimiters;
    
    pos = string::npos;
    linenum = 0;
    
  }

  ~FileParser() {
    
    ifile.close();
  
  }
  
  void skipLine() {
    
    if (pos != string::npos)
      pos = string::npos;
    else {
      getline(ifile,line); 
      ++linenum;
    }
      
  }
    
  FileParserTokenType getNextToken(string& token,const FileParserTokenType expectedTokenType) {
    
    // error checking version...
    
    if (getNextToken(token) != expectedTokenType) {
      cerr << "Error: FileParser type mismatch: token=\"" << token << "\" expectedTokenType=" << expectedTokenType << " line=" << linenum << endl;
      throw(-1);
    }
    
    return(expectedTokenType);
    
  }

  FileParserTokenType getNextToken(string& token) {
    
    // if we came into this call with pos not as npos, then it means there may be
    // more tokens left in the line...
    
    if (pos != string::npos)
      pos = line.find_first_not_of(whitespace,pos);

    // loop until we get something non-whitespace... 

    while (pos == string::npos) {
      getline(ifile,line); 
      ++linenum;
      if (!ifile.good())
	return(FILEPARSER_NO_TOKEN);
      pos = line.find_first_not_of(whitespace,0);
    }
    
    // we have a non-whitespace pos in string line...
    
    //cout << "line: \"" << line << "\" with pos=" << pos << " starting char=\"" << line[pos] << "\"" << endl;
    
    // check if the char is a delimiter...
    
    if (delimiters.find(line[pos]) != string::npos) {

      // if the next char is part of the delimiters, then return
      // just that char in token and let the calling process know...

      token = line[pos++];
      return(FILEPARSER_DELIMITER_TOKEN);

    }
    else {

      // otherwise it is a standart token and 
      
      const string::size_type pos0 = pos;
      pos = line.find_first_of(whitespace_and_delimiters,pos0);
      token = line.substr(pos0,pos-pos0);
      return(FILEPARSER_STANDARD_TOKEN);
      
    }
    
  }

  int stringToInt(const string& token) const {
    int value;
    if (!from_string<int>(value,token,std::dec)) {
      cerr << "Error: FileParser cannot convert string to int: \"" << token << "\"" << endl;
      throw(-1);
    }
    return(value);
  }

  double stringToDouble(const string& token) const {
    double value;
    if (!from_string<double>(value,token,std::dec)) {
      cerr << "Error: FileParser cannot convert string to double: \"" << token << "\"" << endl;
      throw(-1);
    }
    return(value);
  }

  

    
  
};

#endif
