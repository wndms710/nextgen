#ifndef UGRID_READER_HPP
#define UGRID_READER_HPP

#include <math.h>
#include <assert.h>
#include "ByteSwap.hpp"
#include "AsciiReader.hpp"

using namespace std;

// Information on the file format is taken from:
// http://www.simcenter.msstate.edu/software/downloads/doc/ug_io/3d_grid_file_type_ugrid.html

class UgridReader {
protected:
  int counts[7];       // header counts
  int n_zones;

  UgridReader() {
    for (int i=0; i<7; ++i) counts[i] = 0;
    n_zones = 0;
  }

  ~UgridReader() {}

public:
  // methods to be implemented by derived classes if needed
  virtual void readHeader(int& nsp,int& nst) {};
  virtual void readNodes(double (*xsp)[3]) {};
  virtual void readTrisAndQuads(int (*spost)[3],int* znost) {};

  // used by binary
  virtual void setFortran() {};
  virtual void setBytes(const int) {};
  virtual void setByteSwap() {};
};

#define UGRID_FORTRAN_OFFSET  4
#define UGRID_INT_BYTES       4

class UgridBinaryReader : public UgridReader {
private:

  int n_bytes;         // number of bytes per float/double value in binary file
  long int meshOffsets[9];  // file offset for different data sections
  bool byteSwap;       // big vs. little endian
  bool b_fortran;

protected:
  FILE * fp;

public:

  UgridBinaryReader() {
    byteSwap = false;
    b_fortran = false;
    fp = NULL;
  }

  UgridBinaryReader(const string& filename) : UgridReader() {
    init(filename);
  }

  ~UgridBinaryReader() {
    if (fp != NULL) fclose(fp);
  }

  void init(const string& filename) {
    n_bytes = 4;
    byteSwap = false;
    b_fortran = false;

    const int file_err = MiscUtils::openFile(&fp,filename,"rb");
    if (file_err != 0) throw(file_err);
  }

  void readHeader(int& nsp,int& nst) {

    rewindFile();  // start at beginning of file

    // fortran includes 4-byte offset before and after sectioned writes
    if (b_fortran) fseek (fp,UGRID_FORTRAN_OFFSET,SEEK_CUR);
    fread(&counts,UGRID_INT_BYTES,7,fp);
    if (b_fortran) fseek (fp,UGRID_FORTRAN_OFFSET,SEEK_CUR);

    if (byteSwap) {
      for (int i=0; i<7; ++i) counts[i] = ByteSwap::byteSwap(counts[i]);
    }
    COUT1(" > header counts:\n"
      << "    > nodes             : " << counts[0] << "\n"
      << "    > surface triangles : " << counts[1] << "\n"
      << "    > surface quads     : " << counts[2] << "\n"
      << "    > 4-node tet   (vol): " << counts[3] << "\n"
      << "    > 5-node tet   (vol): " << counts[4] << "\n"
      << "    > 6-node tet   (vol): " << counts[5] << "\n"
      << "    > hex elements (vol): " << counts[6]);

    nsp = counts[0];
    nst = counts[1] + 2*counts[2];  // tris then quads (we split these ourselves)

    // compute file offsets:

    // start of node section
    meshOffsets[0] = UGRID_INT_BYTES*7;
    if (b_fortran) {
      meshOffsets[0] += 3*UGRID_FORTRAN_OFFSET;  // account for header open & close, and mesh data section open
    }

    // start of surface tris
    meshOffsets[1] = meshOffsets[0] + counts[0]*n_bytes*3;  // node record

    // start of surface quads
    meshOffsets[2] = meshOffsets[1] + counts[1]*UGRID_INT_BYTES*3;  // tri nodes

    // start of surface zone record
    meshOffsets[3] = meshOffsets[2] + counts[2]*UGRID_INT_BYTES*4;  // quad nodes

    // start of volume 4-nodes
    meshOffsets[4] = meshOffsets[3] + (counts[1]+counts[2])*UGRID_INT_BYTES;  // boundary zones

    // start of volume 5-nodes
    meshOffsets[5] = meshOffsets[4] + (counts[3])*UGRID_INT_BYTES*4;  // 4-nodes

    // start of volume 6-nodes
    meshOffsets[6] = meshOffsets[5] + (counts[4])*UGRID_INT_BYTES*5;  // 5-nodes

    // start of volume hex
    meshOffsets[7] = meshOffsets[6] + (counts[5])*UGRID_INT_BYTES*6;  // 6-nodes

    // start of optional section
    meshOffsets[8] = meshOffsets[7] + (counts[6])*UGRID_INT_BYTES*7;  // 8-nodes (hex)
    if (b_fortran) {
      meshOffsets[8] += 1*UGRID_FORTRAN_OFFSET;  // account for mesh data section close
    }
  }

  void readNodes(double (*xsp)[3]) {

    // start after header
    fseek(fp,meshOffsets[0],SEEK_SET);

    double vals[3];
    for (int isp=0; isp<counts[0]; ++isp) {
      fread(vals,n_bytes,3,fp);
      FOR_I3 xsp[isp][i] = (byteSwap) ? ByteSwap::byteSwap(vals[i]) : vals[i];
    }
  }

  void readTrisAndQuads(int (*spost)[3],int* znost) {

    fseek(fp,meshOffsets[1],SEEK_SET);

    // first triangle spost
    {
      int vals[3];
      for (int ist=0; ist<counts[1]; ++ist) {
        fread(vals,UGRID_INT_BYTES,3,fp);
        FOR_I3 spost[ist][i] = (byteSwap) ? ByteSwap::byteSwap(vals[i]) : vals[i];
      }
    }

    // now quad spost
    // we will split these ourselves
    {
      int vals[4];
      for (int ist=0,end=2*counts[2]; ist<end; ist+=2) {
        fread(vals,UGRID_INT_BYTES,4,fp);

        // first tri
        spost[counts[1]+ist][0] = (byteSwap) ? ByteSwap::byteSwap(vals[0]) : vals[0];
        spost[counts[1]+ist][1] = (byteSwap) ? ByteSwap::byteSwap(vals[1]) : vals[1];
        spost[counts[1]+ist][2] = (byteSwap) ? ByteSwap::byteSwap(vals[2]) : vals[2];

        // second tri
        spost[counts[1]+ist+1][0] = (byteSwap) ? ByteSwap::byteSwap(vals[0]) : vals[0];
        spost[counts[1]+ist+1][1] = (byteSwap) ? ByteSwap::byteSwap(vals[2]) : vals[2];
        spost[counts[1]+ist+1][2] = (byteSwap) ? ByteSwap::byteSwap(vals[3]) : vals[3];
      }
    }

    // now read zone info
    // triangles
    set<int> unique_zones;  // count
    int vals;
    for (int ist=0; ist<counts[1]; ++ist) {
      fread(&vals,UGRID_INT_BYTES,1,fp);
      znost[ist] = (byteSwap) ? ByteSwap::byteSwap(vals) : vals;
      unique_zones.insert(znost[ist]);
    }
    // quads
    for (int ist=0,end=2*counts[2]; ist<end; ist+=2) {
      fread(&vals,UGRID_INT_BYTES,1,fp);

      // first tri
      znost[counts[1]+ist] = znost[counts[1]+ist+1] = (byteSwap) ? ByteSwap::byteSwap(vals) : vals;
      unique_zones.insert(znost[counts[1]+ist]);
    }

    COUT1(" > number of zones described in file: " << unique_zones.size());
  }

  void setFortran() {
    b_fortran = true;
  }

  void setBytes(const int n) {
    assert(n==4 || n==8);
    n_bytes = n;
  }

  void setByteSwap() {
    byteSwap = true;
  }

  void rewindFile() {
    rewind(fp);
  }

};

#undef UGRID_FORTRAN_OFFSET
#undef UGRID_INT_BYTES

enum ugridAsciiPos {
  START,
  AFTER_HEADER,
  AFTER_NODES,
  AFTER_TRIS,
  AFTER_QUADS,
  AFTER_ZONES,
  AFTER_TETS,
  AFTER_PENT5,
  AFTER_PENT6,
  AFTER_HEX
};

class UgridAsciiReader : public UgridReader, public AsciiReader {
private:
  ugridAsciiPos file_pos;

public:

  UgridAsciiReader() {
    file_pos = START;
  }

  UgridAsciiReader(const string& filename) : UgridReader(), AsciiReader(filename) {
    file_pos = START;
  }

  ~UgridAsciiReader() {}

  void readHeader(int& nsp,int& nst) {
    rewindFile();  // start at beginning of file

    assert(file_pos == START);

    for (int i=0; i<7; ++i) {
      counts[i] = getNextTokenAsInt();
    }
    file_pos = AFTER_HEADER;

    COUT1(" > header counts:\n"
      << "    > nodes             : " << counts[0] << "\n"
      << "    > surface triangles : " << counts[1] << "\n"
      << "    > surface quads     : " << counts[2] << "\n"
      << "    > 4-node tet   (vol): " << counts[3] << "\n"
      << "    > 5-node tet   (vol): " << counts[4] << "\n"
      << "    > 6-node tet   (vol): " << counts[5] << "\n"
      << "    > hex elements (vol): " << counts[6]);

    nsp = counts[0];
    nst = counts[1] + 2*counts[2];  // tris then quads (we split these ourselves)
  };

  void readNodes(double (*xsp)[3]) {
    assert(file_pos = AFTER_HEADER);
    for (int isp=0; isp<counts[0]; ++isp) {
      FOR_I3 xsp[isp][i] = getNextTokenAsDouble();
    }
    file_pos = AFTER_NODES;
  };

  void readTrisAndQuads(int (*spost)[3],int* znost) {
    assert(file_pos = AFTER_NODES);

    // tris
    for (int ist=0; ist<counts[1]; ++ist) {
      FOR_I3 spost[ist][i] = getNextTokenAsInt();
    }

    // quads
    for (int ist=0,end=2*counts[2]; ist<end; ist+=2) {
      int nodes[4];
      for (int i=0; i<4; ++i) nodes[i] = getNextTokenAsInt();
      // first tri
      spost[counts[1]+ist][0] = nodes[0];
      spost[counts[1]+ist][1] = nodes[1];
      spost[counts[1]+ist][2] = nodes[2];

      // second tri
      spost[counts[1]+ist+1][0] = nodes[0];
      spost[counts[1]+ist+1][1] = nodes[2];
      spost[counts[1]+ist+1][2] = nodes[3];
    }

    // zone info
    set<int> unique_zones;  // count
    for (int ist=0; ist<counts[1]; ++ist) {
      znost[ist] = getNextTokenAsInt();
      unique_zones.insert(znost[ist]);
    }
    for (int ist=0,end=2*counts[2]; ist<end; ist+=2) {
      znost[counts[1]+ist] = znost[counts[1]+ist+1] = getNextTokenAsInt();
      unique_zones.insert(counts[1]+ist);
    }

    file_pos = AFTER_ZONES;

    COUT1(" > number of zones described in file: " << unique_zones.size());
  };

  void rewindFile() {
    AsciiReader::rewindFile();
    file_pos = START;
  }
};

class BcAsciiReader : public AsciiReader {

public:

  // boundary condition file specific parsing
  void getZonenameMap(map<int,string>& ugridzone2name,const string filetype) {
    if (filetype == "SLUGG") {
      getSluggZonenameMap(ugridzone2name);
    }
    else {
      WUI(WARN,"only SLUGG file formats are currently supported");
    }
  }

  void getSluggZonenameMap(map<int,string>& ugridzone2name) {
    rewindFile();

    goToNextLine(true); // knowing we are at the start of the line, force moving to the next one
    goToNextLine(true);

    while (true) {

      string zonename = getNextTokenAsString();
      if (b_eof) break;

      zonename = zonename.substr(1,string::npos);  // skip "#"
      goToNextLine();  // skipping other items
      if (b_eof) break;

      const int ugridzone = getNextTokenAsInt();

      ugridzone2name.insert(pair<int,string> (ugridzone,zonename));
      COUT2("    > ugrid zone " << ugridzone << " : " << zonename);
      goToNextLine();  // skipping other items

      if (b_eof) break;
    }
  }

};

#endif
