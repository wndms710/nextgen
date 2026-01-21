#ifndef RESTARTHASHUTILITIES_HPP
#define RESTARTHASHUTILITIES_HPP

#include "Common.hpp"
#include "tomcrypt.hpp"
#include "MpiStuff.hpp"
#include "Params.hpp"
using namespace MpiStuff;

namespace RestartHashUtilities {

  const int sha1hashlength = 20;
 
  int writeHashToRestartHeader(MPI_File &fh);
  void mlesReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header);
  void slesReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header);
  int slesConsistencyCheck();
  void slesSetHashFromParam(Param * param);
  void mlesInterpReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header);
  void slesInterpReadHashes(MPI_File &fh, const MPI_Offset &offset, const Header &header);
  int slesInterpConsistencyCheck();
  int slesConsistencyCheck();

  void readHashFromRestartHeader_s(FILE * fp, const Header &header);
  void mlesReadHashes_s(FILE * fp, const Header &header);
  void slesReadHashes_s(FILE * fp, const Header &header);

  void clearSlesHash();
  void clearHashes();

  //void generateEmptyRestartHash(unsigned char * out);
  //bool isEmptyRestartHash(const unsigned char * const in);


  class RestartHash {
   
    private:
  
      unsigned char * binaryHash;
      int length;
      hash_state sha;
  
    public:
  
      RestartHash();
      ~RestartHash();

      //read and init
      void init(MPI_File &fh, MPI_Offset &offset, const int &_length);  
      void init(FILE * fp, const int &_length);  
      //generate new hash
      void init(const stringstream &_ss,const int &_length);
 
      void bCast0ToAll(); //copy binaryHash from rank0 to all ranks

      void clear();
      unsigned char * getBinaryHash() const;
      int getLength() const;
   
      string getAsciiHash();

      bool operator==(const RestartHash &other) const;
  
    private:
      void initSha1();
      void processHashInput(const char * in, const int &n);
      void generateHash();
  
  };

  void bcastHash0ToAll(RestartHash &rHash);

  extern RestartHash mlesHash;
  extern RestartHash slesHash;
  extern RestartHash myHash;

  extern RestartHash mlesInterpHash;
  extern RestartHash slesInterpHash;

  ostream& operator<<(ostream& out, RestartHash &myRestartHash);
}

#endif
