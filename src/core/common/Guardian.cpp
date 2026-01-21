// ======================================================================
// Guardian.cpp
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

#include "Params.hpp"
#include "Guardian.hpp"
#include <pwd.h>
#include <unistd.h> 


using namespace std;

const string Guardian::ctiLicenseFile = "./license.dat";

#ifndef NO_GUARDIAN

Guardian::Guardian(){
  myLicError = 0; //set to 1 if license checking fails
  licError   = 0; //Collective sum of errors across ranks

  license_data = NULL;
  signature_data = NULL;

  //Import Cascade Public Key
  isKeySet = false;
  if (importKey()==EXIT_FAILURE){
    myLicError=1;
  }
  printLicenseError("There was a problem during the license key import");
  
  //Get System Fingerprint
  if (lookupSystemFingerprint()==EXIT_FAILURE){
    myLicError = 1;
  }
  printLicenseError("Unable to generate system fingerprint");
}

Guardian::~Guardian(){
  DELETE(license_data);
  DELETE(signature_data);
  if (isKeySet)
    rsa_free(&rsaKey);
}


void Guardian::getUName(string& username){
    try {
      struct passwd *pw;
      uid_t uid;
      uid = geteuid();
      pw = getpwuid(uid);
      if (pw) username = pw->pw_name;
      else username = "N/A";
    }
    catch (...){
      username = "N/A";
    }
    sysUN = username;
}

void Guardian::getSerialNumber(string& serialNo) {
    try {
      serialNo = "00000000-0000-0000-0000-000000000000";
    }
    catch (...){
      serialNo = "N/A";
    }
}

void Guardian::initialize(string &filePath){

  //Read in license.dat, throw error if unable
  //set signature_data, license_data
  license_data = NULL;
  signature_data = NULL;
  readLicenseFile(filePath);
  if (license_data == NULL){
    myLicError = 1;
  }
  printLicenseError("Data not found in license file");
  if (signature_data == NULL){
    myLicError = 1;
  }
  printLicenseError("Signature not found in license file");


  //Use public key to check license is authentic, if not throw an error
  //EXIT_FAILURE is returned only if there is a code problem,
  //license problems throw an error
  if (validateLicense()==EXIT_FAILURE){
    myLicError = 1;
  }
  printLicenseError("There was a problem during license validation");

  setLicenseInfo();
}


bool Guardian::isValid(){
  bool valid = true;
  if(!checkExpiration()){
    myLicError = 1;
    valid = false;
  }
  string licExpMsg = "The license expired on " + getLicExp();
  printLicenseError(licExpMsg );

  if(!checkFingerprint()){
    myLicError = 1;
    valid = false;
  }
  printLicenseError("Warning: Failed System ID Check");


  return valid;
}

//bool Guardian::isValid(string &solverName){
//  bool valid = true;
//  if (!checkSolvers(solverName)){
//    myLicError = 1;
//    valid = false;
//  }
//  printLicenseError("This solver is not supported by the current license");
//  valid = isValid();
//
//  return valid;
//}

int Guardian::importKey(){

  unsigned char * public_key = (unsigned char*) "MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAyNh0EAieQuoLSYd5Y9BK\nbUYkXJoKvBfqSiPK1dLsh+wrqLh3YtWkq5ZYtmdcxkcqwpIYaixqp4Ro4haVPBqY\nBtbvU7wMhJYqzNFC+vyB1Hlg2VMcKvAIukjZHkmaIjaFjcC+Kf1qnMxIJJKdgHxi\nkTN1zsmpyPO3zIlmygibUAFfxcrwqf1QjtNAqZ/94OFFs0UVZp29hXljSifB9eLd\nnsZ5guFGAsSNF7P8a8NMWF02uFtEo5qsoBYTx09PEJYdDcPtE6siq/ERO2/Rne1h\n/GuiEr52Eq9RyKcw8SdJQVrjSD6RMDnKaWwdy3rRhkNypVlJ2/hpG4+KPQG3FBOX\nQwIDAQAB";

  /* register a math library (in this case TomsFastMath) */
  ltc_mp = ltm_desc;

  /* 2. import RSA public key - This is the public key. Private key must be in DER form and base64-encoded.*/
  unsigned long public_key_length = strlen((char*)public_key);

  /* Base64 decode to recover DER form - The key is 2048 length */
  unsigned char decoded_key[2048];
  unsigned long decoded_key_length = 2048;
 
  int err;
  if ((err = base64_decode(public_key, public_key_length, decoded_key, &decoded_key_length)) != CRYPT_OK) {
    //printf("base64_decode %s", error_to_string(err)); 
    return EXIT_FAILURE;
  }

  if ((err = rsa_import(decoded_key, decoded_key_length, &rsaKey)) != CRYPT_OK) {
    //printf("base64_decode %s", error_to_string(err)); 
    return EXIT_FAILURE;
  }
  
  //Success
  isKeySet = true;
  return 0;
}

int Guardian::lookupSystemFingerprint(){

  // generate fingerprint from current username and code serial number
  string username, serialNo;

  foundUsername = getFingerprintUNAndSN(username,serialNo);

  char fingerprint[256];
  sprintf(fingerprint, "%s:%s", username.c_str(), serialNo.c_str());

  //Begin hashed fingerprint
  unsigned char hash_msg[21]; // sha1 is 20 char
  hash_state md;
  sha1_init(&md);
  sha1_process(&md, (unsigned char*) fingerprint, (long)strlen(fingerprint));
  sha1_done(&md, hash_msg);
  hash_msg[20] = '\0'; // strip out trash character
  // Base64 encode the digest
  unsigned long encoded_fnghash_length = 2048;
  unsigned char encoded_fnghash[2048+1];
  int err;
  if ((err = base64_encode((unsigned char*)hash_msg, 20, encoded_fnghash, &encoded_fnghash_length)) != CRYPT_OK) {
    printf("base64_decode %s", error_to_string(err)); return EXIT_FAILURE;
  }
  encoded_fnghash[encoded_fnghash_length] = '\0';

  sysFngHash = (char*) encoded_fnghash;
  //End hashed fingerprint
  
  //Begin encoded fingerprint
  ltc_mp = ltm_desc;  //math setup
  int hash_idx, prng_idx;
  hash_idx = 0; //sha1 by default
  prng_idx = 0; //sprng by default
  memcpy(&hash_descriptor[hash_idx], &sha1_desc, sizeof(struct ltc_hash_descriptor));
  memcpy(&prng_descriptor[prng_idx], &sprng_desc, sizeof(struct ltc_prng_descriptor));

  //Encrypt the fingerprint with Public Key
  unsigned char encrypted_text[4096];
  unsigned long encrypted_text_length = 4096;
  if ((err = rsa_encrypt_key_ex((unsigned char*)fingerprint, strlen(fingerprint), encrypted_text, &encrypted_text_length, NULL, 0, NULL, prng_idx, hash_idx, 1, &rsaKey)) != CRYPT_OK) {
    return EXIT_FAILURE;
  }

  // Base64 encode the digest 
  unsigned long encoded_fingerprint_length = 2048;
  //unsigned char encoded_fingerprint[encoded_fingerprint_length+1];
  unsigned char encoded_fingerprint[2048+1];

  if ((err = base64_encode((unsigned char*)encrypted_text, encrypted_text_length, encoded_fingerprint, &encoded_fingerprint_length)) != CRYPT_OK) {
    return EXIT_FAILURE;
  }
  encoded_fingerprint[encoded_fingerprint_length] = '\0';

  //END encoded fingerprint

  sysFng = (char *) encoded_fingerprint;

  //returning zero indicates success, EXIT_FAILURE == 1
  return 0;
}

void Guardian::readLicenseFile(const string& filename){
  
  //Delete if arrays exist, set to NULL
  DELETE(license_data);
  DELETE(signature_data);
  
  myLicError = 1;
  string line;
  ifstream myfile(filename.c_str());
  if (myfile.is_open()){
    while ( getline(myfile,line) ){
      string token = line.substr(0,12);
      //cout << " TOKEN: " << token << endl;
      if(token=="#CTI_LICENSE"){
        getline(myfile,line);
        int line_length = strlen(line.c_str());
        //License data is first line of file
        assert(license_data == NULL);
        license_data = new unsigned char[line_length+1];
        memcpy(license_data, line.c_str(), line_length+1);
        //Signature data is the second line of file
        getline(myfile,line);
        line_length = strlen(line.c_str());
        assert(signature_data == NULL);
        signature_data = new unsigned char[line_length+1];
        memcpy(signature_data, line.c_str(), line_length+1);
        myLicError = 0;
        break;
      }
    }
    myfile.close();
  }
  printLicenseError("License file not found, check CTI_LICENSE");
}

int Guardian::validateLicense(){

  int err, hash_idx, prng_idx, res;

  /* register a math library (in this case TomsFastMath) */
  ltc_mp = ltm_desc;

  /* register hash and prng */
  hash_idx = 0; //sha1 by default
  prng_idx = 0; //sprng by default
  memcpy(&hash_descriptor[hash_idx], &sha1_desc, sizeof(struct ltc_hash_descriptor));
  memcpy(&prng_descriptor[prng_idx], &sprng_desc, sizeof(struct ltc_prng_descriptor));

  /* All data is encoded to Base64 for transferring */

  /* 1. Base64 decode for signature data */
  unsigned long signature_data_length = strlen((char*)signature_data);

  unsigned char decoded_signature_data[1024];
  unsigned long decoded_signature_data_length = 1024;
  if ((err = base64_decode(signature_data, signature_data_length, decoded_signature_data, &decoded_signature_data_length)) != CRYPT_OK) {
      myLicError = 1;
//    printf("base64_decode %s", error_to_string(err)); 
//    return EXIT_FAILURE;
  }
  printLicenseError("Unable to decode license signature");

  /* 3. Verify the signature data */

  /* Create a SHA1 hash from license data first */
  unsigned long license_data_length = strlen((char*)license_data);

  unsigned char hash_msg[21]; //
  hash_state md;
  sha1_init(&md);
  sha1_process(&md, license_data, license_data_length);
  sha1_done(&md, hash_msg);
  hash_msg[20] = '\0'; // strip out trash character

  unsigned long hash_msg_length = 20;

  /* Verify hash with signature */
  if ((err = rsa_verify_hash_ex(decoded_signature_data, decoded_signature_data_length, hash_msg, hash_msg_length, 1, hash_idx, 0, &res, &rsaKey)) != CRYPT_OK) {
    myLicError = 1;
//    printf("base64_decode %s", error_to_string(err));
//    return EXIT_FAILURE;
  }
  printLicenseError("Unable to verify license data hash with signature");


  if (!(res == 1)){
    myLicError = 1;
//    return EXIT_FAILURE;
  }
  printLicenseError("Non-matching license signature");



  //Decode License data and store JSON to string
  //JSON string may not be longer than 2047 characters
  const int jsonMaxLength = 2048;
  unsigned char decoded_lic[jsonMaxLength];
  unsigned long decoded_lic_length = jsonMaxLength;
  int err_temp;
  if ((err_temp = base64_decode(license_data, strlen((char*) license_data), decoded_lic, &decoded_lic_length)) != CRYPT_OK) {
    myLicError = 1;
//    printf("base64_decode %s", error_to_string(err_temp)); 
//    return EXIT_FAILURE;
  }
  printLicenseError("Unable to decode license data");

  if (decoded_lic_length < jsonMaxLength)
    decoded_lic[decoded_lic_length] = '\0';  
  else{
    myLicError = 1;
//    return EXIT_FAILURE;
  }
  printLicenseError("License data too long for system buffer");
    
  licJson = (char*) decoded_lic;
  //cout << "DAP LicJson: \n " << licJson << endl;

  //returning zero indicates success, EXIT_FAILURE == 1
  return 0;
}

void Guardian::setLicenseInfo(){

  //Expected license fields, parse licJson string
  map<string,string> licMap;
  licMap["organization"] = "str";
  licMap["username"]     = "str";
  licMap["serial"]       = "str";
  licMap["expiration"]   = "int";
  licMap["fingerprint"]  = "str";

  for (map<string,string>::iterator it=licMap.begin(); it!=licMap.end(); ++it){
    stringstream ss_key;
    ss_key << "\"" << it->first << "\"";
    int i_key = licJson.find(ss_key.str());
    int n_key = ss_key.str().length();
    
    if (i_key!=string::npos){
      int i0 = 0;
      int i1 = 0;
      if (it->second=="str"){
        i0 = licJson.find_first_of("\"",i_key+n_key);
        if (i0!=string::npos)
          i1 = licJson.find_first_of("\"",i0+1);
      }
      else if (it->second=="int") {
        i0 = licJson.find_first_of(":",i_key+n_key);
        if (i0!=string::npos)
          i1 = licJson.find_first_of(",}",i0);
      }
      else {
        myLicError = 1;
        printLicenseError("Unexpected JSON data type");
      }

      if (i0!=string::npos&&i1!=string::npos){
        it->second = licJson.substr(i0+1,i1-i0-1);
      }
      else{
        it->second = ""; //problem parsing, blank this field
        myLicError = 1;
        printLicenseError("Problem parsing license JSON data");
      }
    }
    else{
      it->second = ""; //if key is not found, blank this field
      //Permit fields to be missing, failure will occur later
      //if it is not simply used for reporting
      //myLicError = 1;
      //string errorMsg = "Missing license JSON key: " + it->first;
      //printLicenseError(errorMsg);
    }
    //cout << "licMap[" << it->first << "] = " << it->second << endl;
  }

  licUsr     = licMap["username"];
  licOrg     = licMap["organization"];
  licFngHash = licMap["fingerprint"];
  licSer     = licMap["serial"];
  stringstream exp_ss(licMap["expiration"]);
  exp_ss >> licExp;

//  // Should be an array
//  if ( !licValue.isJBoxObject() ){
//    myLicError = 1;
//  }
//  printLicenseError("Invalid JSON syntax in license data");
//
}

bool Guardian::checkExpiration(){
  time_t now = time(0);
  //Adjust license time for user's timezone
  struct tm * timeinfo;
  timeinfo = gmtime(&licExp);
  time_t localExp = mktime(timeinfo);
  return (now <= localExp);
}

bool Guardian::checkFingerprint(){
#ifdef WITH_RESCALE
  foundUsername = false;
#endif
  return ( !foundUsername || (licFngHash == sysFngHash) );
}

//bool Guardian::checkSolvers(string &solverKey){
//  //All license files should have the
//  //'ALL' entry as true or false
//  bool isValidSolver = licSolvers["ALL"].getBoolean();
// 
//  if (!isValidSolver){
//    JBoxObject::iterator objIt = licSolvers.find(solverKey);
//    if (objIt != licSolvers.end()){
//       isValidSolver = objIt->second.getBoolean();
//    }    
//  }
//
//  return isValidSolver;
//}

//string Guardian::listAllowedSolvers(){
//  string theList = "";
//  int i = 0;
//  if (licSolvers["ALL"].getBoolean())
//    theList += "All";
//  else {
//    for (JBoxObject::iterator objIt = licSolvers.begin(); objIt != licSolvers.end(); ++objIt)
//      if (objIt->second.getBoolean()){
//        if (i>0){
//          theList += ", ";
//        }
//        theList += objIt->first;
//        i++;
//      }
//  }
//  return theList;
//}

string Guardian::getLicUsr() const {
  return licUsr;
}

string Guardian::getLicOrg() const {
  return licOrg;
} 

string Guardian::getLicExp() const {
  //Adjust license time for user's timezone
  struct tm * timeinfo;
  timeinfo = gmtime(&licExp);
  time_t localExp = mktime(timeinfo);
  return ctime(&localExp);
}

string Guardian::getLicSer() const {
  return licSer;
}

string Guardian::getSysFng() const {
  return sysFng;
}

//Virtual Function serial implementations
  bool Guardian::getFingerprintUNAndSN(string &username,string &serialNo){

    getUName(username);
    getSerialNumber(serialNo);

    //If unable to get username,
    //skip system fingerprint check
    foundUsername=true;
    if (username == "N/A"){
      foundUsername = false;
      licError = 0;
      myLicError = 0;
    }
    return foundUsername;
  }

  void Guardian::printLicenseError(const string &message){
    if (myLicError>0){ 
        if (foundUsername){
        std::cout <<                   
        "\n\n************************ LICENSE ERROR *****************************\n > " << 
        message << "\n" << 
        "\n       The following resource username can be used to generate" << 
        "\n           a license at support.cascadetechnologies.com\n\n" << 
          "             > " << sysUN <<
        "\n\n********************* END OF LICENSE ERROR *************************\n" << std::endl; 
        }
        else{
        std::cout <<                   
        "\n\n************************ LICENSE ERROR *****************************\n > " << 
        message << 
        "\n********************* END OF LICENSE ERROR *************************\n" << std::endl;
        }
      throw(0);
    }
  }

  void Guardian::initialize(){
    //Set path for license file, order preference:
    // env var, compilation default
    bool foundPath = false;
    string filePath;
    if (!foundPath){
      char* tPath;
      tPath = getenv("CTI_LICENSE");
      if (tPath!=NULL){
        filePath = tPath;
        cout << " > Found env CTI_LICENSE=" << filePath << endl;;
        foundPath = true;
      }
    } 

    if (!foundPath){
      filePath = ctiLicenseFile;
      cout << " > Using default path, CTI_LICENSE=" << filePath << endl;
    }

    Guardian::initialize(filePath);
  }

  bool Guardian::isFingerprintCheckEnabled() const{ 
    return foundUsername;
  }

#else
 
Guardian::Guardian(){
  myLicError     = 0; 
  licError       = 0; 
  license_data   = NULL;
  signature_data = NULL;
  isKeySet       = false;
}

Guardian::~Guardian(){
  DELETE(license_data);
  DELETE(signature_data);
}


void Guardian::getUName(string& username){} 
void Guardian::getSerialNumber(string& serialNo) {} 

void Guardian::initialize(string &filePath){} 

bool Guardian::isValid(){ return true; } 
//bool Guardian::isValid(string &solverName){ return true;} 

int Guardian::importKey(){ return 0; } 
int Guardian::lookupSystemFingerprint(){ return 0;} 

void Guardian::readLicenseFile(const string& filename){} 
int Guardian::validateLicense(){ return 0;} 
void Guardian::setLicenseInfo(){} 

bool Guardian::checkExpiration(){ return true;} 
bool Guardian::checkFingerprint(){ return true;} 
//bool Guardian::checkSolvers(string &solverKey){ return true;} 
//string Guardian::listAllowedSolvers(){ return "";} 

string Guardian::getLicUsr() const { return licUsr; } 
string Guardian::getLicOrg() const { return licOrg; } 
string Guardian::getLicExp() const { return ctime(&licExp); } 
string Guardian::getLicSer() const { return licSer; } 
string Guardian::getSysFng() const { return sysFng; } 

bool Guardian::getFingerprintUNAndSN(string &username,string &serialNo){ return true;} 
void Guardian::printLicenseError(const string &message){} 

void Guardian::initialize(){} 

bool Guardian::isFingerprintCheckEnabled() const { return false;}

#endif



