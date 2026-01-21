#include "SimpleRestart.hpp"
#include "Defs.hpp"
#include "ByteSwap.hpp"
#include "Utils.hpp"


void SimpleRestart::clear(){
  filetype = "unknown";
  filename_mesh = "";
  filename_data = "";
  ss.clear();
  meshMaxOffset = -1;
  cvD1OffsetMap.clear();
  cvD2OffsetMap.clear();
}

int SimpleRestart::init(Param * param) {
  COUT1("SimpleRestart::init()");

  int iarg = 0;
  int ierr = 0;
  if (param->size() > 0) {
    filename_mesh = param->getString(iarg++);
    cout << "about to read restart file: " << filename_mesh << endl;
    ierr = initMLES();

    if (ierr == 0){
      bool b_foundSLES = false;
      if (iarg < param->size()) {
        filename_data = param->getString(iarg++);
        cout << "about to read vars from solution file: " << filename_data << endl;
        ierr = initSLES();
        if (ierr != 0) {
          CWARN("could not read sles file; only mesh variables will be available"); //TODO WUI INFO?
          ierr = 0;
        }
        else{
          b_foundSLES = true;
        }
      }
      // determine if valid vars for mesh drawing are available
      bool b_requiredVars = false;
      map<const string,int8>::iterator it1;
      map<const string,int8>::iterator it2;
      it1 = cvD2OffsetMap.find("x_vv");
      it2 = cvD2OffsetMap.find("X_VV");
      if (it1 != cvD2OffsetMap.end() || it2 != cvD2OffsetMap.end()) {
        it1 = cvD1OffsetMap.find("r_vv");
        it2 = cvD1OffsetMap.find("R_VV");
        if (it1 != cvD1OffsetMap.end() || it2 != cvD1OffsetMap.end()){
          b_requiredVars = true;
        }
      }

      if (b_requiredVars!=true){
         CWARN("x_vv and/or r_vv not found in " << filename_mesh << ", post cannot analyze this file");
         ierr = 1;
      }
      else{
        cout << "mles detected, data visualization enabled from " << filename_mesh << endl;
        filetype = "mesh";
        if (b_foundSLES) {
          cout << "sles detected, data visualization enabled from " << filename_data << endl;
          filetype = "result";
        }
      }

      if (ierr==0){
        assert(ss.znost);
        WebUI::webUIOutput.ensureImage();
      }
    }  //if ierr==0
  }  //if param->size()>0
  else{
    ierr = 1;
    CWARN("no restart file name specified");
  }

  //if there was an error, make sure nothing is partially initialized;
  if (ierr!=0)
    clear();

  return ierr;  //0 == SUCCESS
}

int SimpleRestart::initMLES() {
  int ierr = initMLES_surface();
  if (ierr != 0){
    CWARN("unable to read surface from " << filename_mesh);         
  }
  else{
    ierr = initMLES_data();
  }
  return ierr;
}
  
int SimpleRestart::initMLES_surface(){
  Param paramSurface("SURF");
  paramSurface.tokens.push_back("RESTART");
  paramSurface.tokens.push_back(filename_mesh);
  return ss.init(&paramSurface);
}
  

//Add cvD1 and cvD2 data to the respective maps
int SimpleRestart::initMLES_data() {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename_mesh,"rb");
  if (file_err != 0) return file_err;

  int byte_swap = 0;
  int itmp[2];
  fread(itmp,sizeof(int),2,fp);
  if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
    ByteSwap::byteSwap(itmp,2);
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      cerr << "Error: file does not start as expected. aborting."<< endl;
      throw(-1);
    }
    cout << "File requires byte swapping." << endl;
    byte_swap = 1;
  }

  const int io_version = itmp[1];
  cout << " > io version: " << io_version << endl;
  assert(io_version >= 4);
  if (io_version != 5){
    CWARN(" >  Warning, an additional R2 variable (r_vv) may need to be added for io_version<5 for compatibility with surfer2");
  }

  bool b_foundHash = false;
  int8 offset = sizeof(int)*2;
  Header header;

  int done = 0;
  while (done != 1) {

    fseek(fp,offset,SEEK_SET);

    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

    cout << "header \"" << header.name << "\" id: " << header.id << endl;

    //Since expr parsing added to surfer, no longer use lowercase for all map keys
    string headerNameKey = header.name;
    //transform(headerNameKey.begin(),headerNameKey.end(), headerNameKey.begin(), ::tolower);

    switch (header.id) {
      case UGP_IO_HASHIDS:
        RestartHashUtilities::clearHashes();
        RestartHashUtilities::mlesReadHashes_s(fp,header);
        cout << " > mles hash: " << RestartHashUtilities::mlesHash << endl;
        b_foundHash = true;
        break;

      case UGP_IO_CV_D1:
        cvD1OffsetMap[headerNameKey] = offset + sizeof(Header);
        break;

      case UGP_IO_CV_D2:
        cvD2OffsetMap[headerNameKey] = offset + sizeof(Header);
        break;

      case UGP_IO_EOF:
        done = 1;
    }

    offset += header.skip;

  }

  //store max offset in this file, sles offsets will be stored in the same map
  //but shifted by this value.
  meshMaxOffset = offset;

  if (!b_foundHash){ //set a default mles hash if one was not found
    RestartHashUtilities::clearHashes();
    std::stringstream ss;
    ss << "Missing mles Hash";
    RestartHashUtilities::mlesHash.init(ss,RestartHashUtilities::sha1hashlength);
    cout << " > setting default mles hash id " <<  RestartHashUtilities::mlesHash << endl;
  }

  fclose(fp);

  cout << " > done read" << endl;

  return 0; //TODO error checking

}

int SimpleRestart::initSLES() {

  //mles must be read first...
  assert(meshMaxOffset>=0);

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename_data,"rb");
  if (file_err != 0) return file_err;

  int byte_swap = 0;
  int itmp[2];
  fread(itmp,sizeof(int),2,fp);
  // sles files's magic number is UGP_IO_MAGIC_NUMBER+1
  if ( (itmp[0] != UGP_IO_MAGIC_NUMBER) && (itmp[0] != UGP_IO_MAGIC_NUMBER+1) ) {
    ByteSwap::byteSwap(itmp,2);
    if ( (itmp[0] != UGP_IO_MAGIC_NUMBER) && (itmp[0] != UGP_IO_MAGIC_NUMBER+1) ) {
      cerr << "Error: file does not start as expected. aborting."<< endl;
      throw(-1);
    }
    cout << "File requires byte swapping." << endl;
    byte_swap = 1;
  }

  const int io_version = itmp[1];
  cout << " > io version: " << io_version << endl;
  if (io_version != 5){
    CWARN(" >  Warning, an additional R2 variable (r_vv) may need to be added for io_version<5 for compatibility with surfer2");
  }

  // initialize the offset...
  bool b_foundHash = false;
  if (io_version<5)
    b_foundHash = true; //skip sles check for older restart files...
  int8 offset = sizeof(int)*2;
  Header header;
  vector<pair<string,int> > zoneNameAndCount;

  int done = 0;
  while (done != 1) {

    fseek(fp,offset,SEEK_SET);

    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

    cout << "header \"" << header.name << "\" id: " << header.id << endl;

    //Since expr parsing added to surfer, no longer use lowercase for all map keys
    string headerNameKey = header.name;
    //transform(headerNameKey.begin(),headerNameKey.end(), headerNameKey.begin(), ::tolower);

    switch (header.id) {

      case UGP_IO_CV_D1:
        cvD1OffsetMap[headerNameKey] = meshMaxOffset + offset + sizeof(Header);
        break;

      case UGP_IO_CV_D2:
        cvD2OffsetMap[headerNameKey] = meshMaxOffset + offset + sizeof(Header);
        break;

      case UGP_IO_HASHIDS:
        if (io_version>=5){
          b_foundHash = true;
          int slesNonMatchFlag;
          //read two hash id's from sles file. First id identifies
          //the current file, store in slesHash.  Second  id
          //identifies the "parent" mles file.  Store in mlesHash
          RestartHashUtilities::slesReadHashes_s(fp, header);
          slesNonMatchFlag = RestartHashUtilities::slesConsistencyCheck(); // 0 match, 1 no match
          cout << " > Found sles hash: " << RestartHashUtilities::slesHash << endl;
          if (slesNonMatchFlag==0){
            cout << " >  with mles hash: " << RestartHashUtilities::mlesHash << endl;
          }
          else{
            cout << " > sles expects mles hash: " << RestartHashUtilities::myHash << endl;
            cout << " >      current mles hash: " << RestartHashUtilities::mlesHash << endl;
            CERR("sles data file is not a match for the existing mles file");
          }
        }
        break;

      case UGP_IO_EOF:
        done = 1;

        break;

    }

    offset += header.skip;
  }
  fclose(fp);

  if (!b_foundHash){ //set a default sles hash if one was not found
    std::stringstream ss;
    ss << "Missing sles Hash";
    RestartHashUtilities::slesHash.init(ss,RestartHashUtilities::sha1hashlength);
    if (mpi_rank==0) cout << " > setting default sles hash id " <<  RestartHashUtilities::slesHash << endl;
  }
  return 0;

}

bool SimpleRestart::isInit() const {
  return !((ss.xsp == NULL)||(ss.spost==NULL)||(ss.nst==0)||(ss.nsp==0));
}


//if returns true, user must cleanup memory
bool SimpleRestart::readResultParams(char * &cbuf, const string& filename) {
  assert(cbuf==NULL);
 
  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;

  int byte_swap = 0;
  int itmp[2];
  fread(itmp,sizeof(int),2,fp);
  // sles files's magic number is UGP_IO_MAGIC_NUMBER+1
  if ( (itmp[0] != UGP_IO_MAGIC_NUMBER) && (itmp[0] != UGP_IO_MAGIC_NUMBER+1) ) {
    ByteSwap::byteSwap(itmp,2);
    if ( (itmp[0] != UGP_IO_MAGIC_NUMBER) && (itmp[0] != UGP_IO_MAGIC_NUMBER+1) ) {
      cerr << "Error: file does not start as expected. aborting."<< endl;
      throw(-1);
    }
    cout << "File requires byte swapping." << endl;
    byte_swap = 1;
  }
  const int io_version = itmp[1];
  cout << " > io version: " << io_version << endl;

  // initialize the offset...
  int8 offset = sizeof(int)*2;
  Header header;

  bool b_foundParams = false;
  int done = 0;
  while (done != 1) {

    fseek(fp,offset,SEEK_SET);

    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

    cout << "header \"" << header.name << "\" id: " << header.id << endl;

    switch (header.id) {

      case UGP_IO_PARAMS:
        cbuf = new char[header.idata[0]+1];
        fread(cbuf,sizeof(char),header.idata[0],fp);
        cbuf[header.idata[0]] = '\0';
        done = 1;
        b_foundParams = true;
        break;

      case UGP_IO_EOF:
        done = 1;
        break;
    }
    offset += header.skip;
  }
  fclose(fp);

  return b_foundParams;

}



