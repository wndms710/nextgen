#include "RestartHashUtilities.hpp"

namespace RestartHashUtilities {

  RestartHash mlesHash;
  RestartHash slesHash;
  RestartHash myHash;

  RestartHash mlesInterpHash;
  RestartHash slesInterpHash;

  //namespace utilities

  int writeHashToRestartHeader(MPI_File &fh){
    if (mpi_rank==0){
      assert(mlesHash.getBinaryHash()&&myHash.getBinaryHash());
      Header header;
      sprintf(header.name, "UGP_IO_HASHIDS");
      header.id = UGP_IO_HASHIDS;
      header.skip = header_size + sizeof(unsigned char)*2*sha1hashlength;
      header.idata[0] = sha1hashlength;
      header.idata[1] = sha1hashlength;
//      MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
//      MPI_File_write_at(fh,offset+header_size,myHash.getBinaryHash(),sha1hashlength,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
//      MPI_File_write_at(fh,offset+header_size+sizeof(unsigned char)*sha1hashlength,mlesHash.getBinaryHash(),sha1hashlength,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      MPI_File_write(fh,myHash.getBinaryHash(),sha1hashlength,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
      MPI_File_write(fh,mlesHash.getBinaryHash(),sha1hashlength,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
    }
    myHash.clear(); //clear myHash in case multiple sles files are written
    return  header_size + sizeof(unsigned char)*2*sha1hashlength;
  }

  //read the child hash and store in mlesHash
  void mlesReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    if (mpi_rank==0){
      MPI_Offset myOffset = offset + header_size; //leave passed in offset unchanged  
      mlesHash.init(fh,myOffset,header.idata[0]);
    }
    mlesHash.bCast0ToAll();
  }

  //use myHash to store mles "parent" hash from sles file.  Check that this
  //matches the existing mles file
  void slesReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    if (mpi_rank==0){
      MPI_Offset myOffset = offset + header_size; //leave passed in offset unchanged  
      slesHash.init(fh,myOffset,header.idata[0]);
      myOffset += header.idata[0];
      myHash.init(fh,myOffset,header.idata[1]);
    }
    slesHash.bCast0ToAll();
  }

  //sometimes no sles file is used to start a run..build one from an INIT_* param
  void slesSetHashFromParam(Param * param){
    //call on all ranks
    slesHash.clear();
    stringstream ss;
    param->dump(ss);
    slesHash.init(ss,sha1hashlength);
  }

  //ensure that the sles parent hash matches the mles hash.  
  //If the hashes do not match, shutdown.
  //expects the sles parent hash is stored in myHash
  //this routine will not be called for legacy mles/sles files without a hash
  //because they have no hash header record
  int slesConsistencyCheck(){
    assert(mpi_rank==0);
    assert(mlesHash.getLength()!=0 && myHash.getLength()!=0);
    if ( !(mlesHash==myHash) ) {// check that hashes match
      return 1;
    }
    myHash.clear();
    return 0;
  }

  /*Function when reading files for interpolation*/

  void mlesInterpReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    if (mpi_rank==0){
      MPI_Offset myOffset = offset + header_size; //leave passed in offset unchanged  
      mlesInterpHash.init(fh,myOffset,header.idata[0]);
    }
    mlesInterpHash.bCast0ToAll();
  }

  void slesInterpReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    if (mpi_rank==0){
      MPI_Offset myOffset = offset + header_size; //leave passed in offset unchanged  
      slesInterpHash.init(fh,myOffset,header.idata[0]);
      myOffset += header.idata[0];
      myHash.init(fh,myOffset,header.idata[1]);
    }
    slesInterpHash.bCast0ToAll();
  }

  int slesInterpConsistencyCheck(){
    assert(mpi_rank==0);
    assert(mlesInterpHash.getLength()!=0 && myHash.getLength()!=0);
    if ( !(mlesInterpHash==myHash) ) {// check that hashes match
      return 1;
    }
    myHash.clear();
    return 0;
  }

  /*End Interp functions*/

  //read child/parent hashes as they are laid out in header [serial]
  void readHashFromRestartHeader_s(FILE * fp, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    myHash.init(fp,header.idata[0]);
    mlesHash.init(fp,header.idata[1]);
  }

  //read the child hash and store in mlesHash [serial version]
  void mlesReadHashes_s(FILE * fp, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength); 
    mlesHash.init(fp,header.idata[0]);
  }

  //use myHash to store mles "parent" hash from sles file.  Check that this
  //matches the existing mles file [serial version]
  void slesReadHashes_s(FILE * fp, const Header &header){
    assert(header.idata[0] <= sha1hashlength && header.idata[1] <= sha1hashlength);
    slesHash.init(fp,header.idata[0]);
    myHash.init(fp,header.idata[1]);
  }

  void clearSlesHash() { 
    slesHash.clear();
  }

  void clearHashes(){
    mlesHash.clear();
    slesHash.clear();
    myHash.clear();
  }

  //RestartHash class methods
  RestartHash::RestartHash() : binaryHash(NULL), length(0) {}
 
  RestartHash::~RestartHash(){
    if (binaryHash){
      delete[] binaryHash;
      binaryHash = NULL;
    }
  }

  void RestartHash::init(FILE * fp, const int &_length) {
    assert(!binaryHash&&length==0);
    length = _length;
    binaryHash = new unsigned char[length+1];
    fread(binaryHash,sizeof(char),length,fp);
    binaryHash[length] = '\0';
  }

  void RestartHash::init(MPI_File &fh, MPI_Offset &offset, const int &_length) {
    assert(!binaryHash&&length==0);
    length = _length;
    binaryHash = new unsigned char[length+1];
    MPI_File_read_at(fh,offset,binaryHash,length,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
    binaryHash[length] = '\0';
  }
  
  void RestartHash::init(const stringstream &_ss,const int &_length) {
    assert(!binaryHash&&length==0);
    length = _length;
    binaryHash = new unsigned char[length+1];
    initSha1();
    const int hashbuffLength = 1024;
    char hashbuff[hashbuffLength+1];
    int totalLength = _ss.str().length();
    int pos = 0;
    int ncpy = 0;
    while (pos < totalLength)
    {
      ncpy = totalLength-pos;
      if (hashbuffLength<ncpy)
        ncpy = hashbuffLength;
  
      _ss.str().copy(hashbuff,ncpy,pos);
      hashbuff[ncpy+1]='\0';
  
      //cout << hashbuff << endl;
  
      processHashInput(hashbuff, ncpy);
      pos += hashbuffLength;
    }
    generateHash();
  }

  void RestartHash::clear(){
    length = 0;
    if (binaryHash){
      delete[] binaryHash;
      binaryHash = NULL;
    }
  }
 
  //send a hash to all ranks, useful for including hash in output that is written
  //from non-rank 0 process (i.e. probes)
  void RestartHash::bCast0ToAll(){
    if (mpi_rank==0) 
      assert(sha1hashlength==length);
    else{
      clear();
      binaryHash = new unsigned char[sha1hashlength+1]; 
    }
    length = sha1hashlength;
    MPI_Bcast(binaryHash,sha1hashlength+1,MPI_UNSIGNED_CHAR,0,mpi_comm);
  }



  unsigned char * RestartHash::getBinaryHash() const {
    return binaryHash;
  }

  int RestartHash::getLength() const {
    return length;
  }

  string RestartHash::getAsciiHash(){
    if (!binaryHash)
      return "NULL";
    unsigned long encoded_hash_length = 2048;
    unsigned char encoded_hash[2048+1];
    base64_encode(binaryHash, length, encoded_hash, &encoded_hash_length);
    string encodedHashString = reinterpret_cast<char*>(encoded_hash);
    return encodedHashString;
  }
  
  //private class methods
  void RestartHash::initSha1(){ //setup hash function, clears old data
    sha1_init(&sha);
  }
  
  void RestartHash::processHashInput(const char * in, const int &n){ //add chars to be hashed (may be repeated)
    unsigned char * in_unsigned = (unsigned char *) in;
    sha1_process(&sha, in_unsigned, n);
  }
  
  void RestartHash::generateHash(){ //get the resulting hash
    sha1_done(&sha, binaryHash);
    binaryHash[length] = '\0';
  }
 
  ostream& operator<<(ostream& out, RestartHash &myRestartHash){
    out << myRestartHash.getAsciiHash();
    return out;
  }

  bool RestartHash::operator==(const RestartHash &other) const {
     if (length!=other.getLength()){
       return false;
     }
     else {
       unsigned char * otherBinaryHash = other.getBinaryHash();
       for(int i=0;i<length;++i){
         if(binaryHash[i]!=otherBinaryHash[i])
           return false;
       }
     }
     return true;
  }

}

