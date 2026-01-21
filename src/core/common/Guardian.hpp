#ifndef CTIGUARDIAN_HPP
#define CTIGUARDIAN_HPP

// ======================================================================
// Guardian.hpp
// 
//
// History:
//  author: David Philips - Oct 2013
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
#include "tomcrypt.hpp"  

using namespace std;

class Guardian {
private:
  //pubilc key data 
  rsa_key rsaKey;
  bool isKeySet;

  //System fingerprint data
  string sysFng;
  string sysFngHash;

  //License file data
  unsigned char * signature_data;
  unsigned char * license_data;
  string licJson;
  string licUsr;
  string licOrg;
  time_t licExp;
  string licSer;
  string licFngHash;

  //JsonBox::JBoxObject licSolvers;
 
  //private methods for license verification
  int importKey();
  int lookupSystemFingerprint();
  void readLicenseFile(const string& filename);
  int validateLicense();
  void setLicenseInfo();

protected:

  static const string ctiLicenseFile;
   
  int licError, myLicError;
  bool   foundUsername;
  string sysUN;

  void getUName(string& username);
  void getSerialNumber(string& serialNo);
  virtual bool getFingerprintUNAndSN(string &username,string &serialNo);

public:
  Guardian();
  virtual ~Guardian();

  virtual void printLicenseError(const string &message);

  virtual void initialize();
  void initialize(string &filePath);

  bool isValid();
  //bool isValid(string &solverName);

  bool checkExpiration();
  bool checkFingerprint();
  //bool checkSolvers(string &solverKey);
 // string listAllowedSolvers();

  //License info printable to log file
  string getLicUsr() const;
  string getLicOrg() const;
  string getLicExp() const;
  string getLicSer() const;
  string getSysFng() const;
  //flag for fingerprint check print statments
  bool isFingerprintCheckEnabled() const; 
};
#endif
