#include "SimpleSurface.hpp"
#include "DoubleVertex.hpp"
// #include "FluentReader.hpp"
#include "parsers/fluent/CasReader.hpp"

int SimpleSurface::addMsh(const string& filename,const int force_byte_swap) {

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  CasReader fr;
  try {
    fr.init(filename);
  }
  catch (...) {
    WUI(WARN,"could not initialize Case File Reader properly for " << filename);
    return -1;
  }

  if (force_byte_swap == 1) fr.setByteSwap(true);

  int tri_count = 0;
  int tri_internal = 0;
  int tri_parent = 0;
  set<int> skipFaceZoneIds;

  COUT1(" > [pass #1] reading face information for potential triangle decimation");
  while (fr.advanceToLevel(1)) {
    const int section_header = fr.getNextTokenAsInt(1);

    if (section_header == 2) {
      // Grid Dimensionality Section
      // format: (2 ND)
      const int nd = fr.getNextTokenAsInt(1);
      COUT2("    > grid dimensionality: " << nd << "-D");

      if (nd != 3) {
        WUI(WARN,"only 3D grids are supported");
        return -2;
      }
    }
    else if (section_header == 4) {
      // machine configuration; will tell you if little or big Endian
      const int endianFlag = fr.getNextTokenAsInt(2);
      if ((endianFlag != 60) && (force_byte_swap != 0)) {
        fr.setByteSwap(true);
        COUT2(" > byte-swapping binary records");
      }
      // for now don't bother with other flags since we don't know what they are....
    }
    else if (section_header == 10) {
      // Node Section
      // format: (10 (id start end type [ND]) (x y [z]))
      // [.] indicate optional parameters, depending on dimension
      const int id    = fr.getNextTokenAsHex(2);
      const int start = fr.getNextTokenAsHex(2);
      const int end   = fr.getNextTokenAsHex(2);
      const int type  = fr.getNextTokenAsHex(2);
      UNUSED(type);

      if (id == 0) {
        assert(start == 1);
        // declaration section, ND is omitted
        // set the node count for checking
        COUT2("    > global node declaration: " << end << " nodes");

        // intialize new surface points data
        nsp += end;

        growNspData(nsp,ss_nsp0);
      }
    }
    else if (section_header == 13) {
      // Face Section
      // format: (13 (id start end bc-type fa-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int bcType = fr.getNextTokenAsHex(2);

      // only process non-header information
      if (id != 0) {
        const int faType   = fr.getNextTokenAsHex(2);

        if ((end-start+1) == 0) {
          COUT2("    ! face block: id: " << id << " has no faces listed in section; skipping");
          skipFaceZoneIds.insert(id);
          continue;
        }
        else {
          COUT2("    > face block: id: " << id << ", nelements: " << end-start+1 << ", type: " << faType);
        }

        if (faType >= 2 && faType <= 4) {
          int tris = (end - start + 1) * (faType - 2);
          int index;
          for (int i = start; i <= end; ++i) {
            for (int j = 0; j < (faType + 2); ++j) {
              index = fr.getNextTokenAsHex(2);
            }
          }

          if (bcType == 2) {
            COUT2("      ! internal face zone detected; skipping");
            tri_internal += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else if (bcType == 31) {
            COUT2("      ! parent hanging-node section; skipping");
            tri_parent += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else {
            COUT2("      - equivalent to " << tris << " tris");
            tri_count += tris;
          }
        }
        else if (faType == 0 || faType == 5) {
          fr.advanceToLevel(2);
          int no_per_fa;
          int tris = 0;
          int cell_index;

          for (int ifa = start; ifa <= end; ++ifa) {
            no_per_fa = fr.getNextTokenAsHex(2);
            tris += no_per_fa - 2;

            for (int j = 0; j < (no_per_fa + 2); ++j) {
              cell_index = fr.getNextTokenAsHex(2);
            }
          }

          if (bcType == 2) {
            COUT2("      ! internal face zone detected; skipping");
            tri_internal += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else if (bcType == 31) {
            COUT2("      ! parent hanging-node section; skipping");
            tri_parent += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else {
            COUT2("      - equivalent to " << tris << " tris");
            tri_count += tris;
          }
        }
      }
    }
    else if ((section_header == 2013) || (section_header == 3013)) {
      // Binary Face Section
      // format: (13 (id start end bc-type fa-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int bcType = fr.getNextTokenAsHex(2);
      const int faType = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // only process non-header information
      if (id != 0) {
        if ((end-start+1) == 0) {
          COUT2("    ! face block: id: " << id << ", no faces listed in section; skipping");
          skipFaceZoneIds.insert(id);
          continue;
        }
        else {
          COUT2("    > face block: id: " << id << ", nelements: " << end-start+1 << ", type: " << faType);
        }

        if (faType >= 2 && faType <= 4) {
          int tris = (end - start + 1) * (faType - 2);

          int x[faType+2];
          for (int ifa = start; ifa <= end; ++ifa) {
            fr.readBinaryFaces(x,faType);
          }

          if (bcType == 2) {
            COUT2("      ! internal face zone detected; skipping");
            tri_internal += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else if (bcType == 31) {
            COUT2("      ! parent hanging-node section; skipping");
            tri_parent += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else {
            COUT2("      - equivalent to " << tris << " tris");
            tri_count += tris;
          }
        }
        else if (faType == 0 || faType == 5) {
          int no_per_fa;
          int tris = 0;

          for (int ifa = start; ifa <= end; ++ifa) {
            fr.readBinaryInt(no_per_fa);
            tris += no_per_fa - 2;

            int x[no_per_fa+2];
            fr.readBinaryFaces(x,no_per_fa);
          }

          if (bcType == 2) {
            COUT2("      ! internal face zone detected; skipping");
            tri_internal += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else if (bcType == 31) {
            COUT2("      ! parent hanging-node section; skipping");
            tri_parent += end-start+1;
            skipFaceZoneIds.insert(id);
          }
          else {
            COUT2("      - equivalent to " << tris << " tris");
            tri_count += tris;
          }
        }
      }
    }
    else if ((section_header == 2010) || (section_header == 3010)) {
      // Need to properly stride across binary sections
      const int id    = fr.getNextTokenAsHex(2);
      const int start = fr.getNextTokenAsHex(2);
      const int end   = fr.getNextTokenAsHex(2);
      const int type  = fr.getNextTokenAsHex(2);
      int noND = 3;
      if (fr.getLevel() == 2) noND  = fr.getNextTokenAsHex(2);

      UNUSED(id);
      UNUSED(type);
      fr.advanceToLevel(2);

      if (section_header == 3010) {
        // read double precision from binary
        double d[3];
        for (int i = start; i <= end; ++i) fr.readBinaryNodes<double>(d,3);
      }
      else {
        assert(section_header == 2010);
        // read single precision from binary
        double d[3];
        for (int i = start; i <= end; ++i) fr.readBinaryNodes<float>(d,3);
      }
      // ------------------------------------------

    }
    else if ((section_header == 2012) || (section_header == 3012)) {
      // Binary Cell Section
      // format: (2012 (id start end type elem-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int type   = fr.getNextTokenAsHex(2);
      const int elemType = fr.getNextTokenAsHex(2);
      UNUSED(id);
      UNUSED(type);

      // only need to read cell description if of mixed-type
      if (elemType == 0) {
        fr.advanceToLevel(2);

        int cellType;
        for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(cellType);
        }
      }
      // otherwise cell block is empty

    }
    else if (section_header == 2040) {
      // Need to properly stride across binary sections
      // (2040 (zone-id first-index last-index partition-count) ())
      const int zone_id  = fr.getNextTokenAsHex(2);
      const int start    = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int part_count    = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if (section_header == 2041) {
      // Need to properly stride across binary sections
      // (2041 (start end) ())
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if (section_header == 2018 || section_header == 3018) {
      // Need to properly stride across binary sections
      // (2041 (start end) ())
      const int start = fr.getNextTokenAsHex(2);
      const int end = fr.getNextTokenAsHex(2);
      const int zone_id  = fr.getNextTokenAsHex(2);
      const int shadow_zone_id = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        // face correspondence between shadow zones
        fr.readBinaryInt(tmp_flag);
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if ((section_header == 2058) || (section_header == 2059)) {
      // Need to properly stride across binary sections
      // (2059/58 (start end parent-zone child-zone) ())
      // where zones are either face or cell blocks, depending
      const int start = fr.getNextTokenAsHex(2);
      const int end = fr.getNextTokenAsHex(2);
      const int parentZone = fr.getNextTokenAsHex(2);
      const int childZone = fr.getNextTokenAsHex(2);
      UNUSED(parentZone);
      UNUSED(childZone);

      fr.advanceToLevel(2);

      int n_children;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(n_children);
        for (int j=0; j<n_children; ++j) fr.getNextBinary<int,int>();  // don't need to store this info
      }
    }
  }
  COUT1("    > triangulation count: " << tri_count);

  fr.rewindFile();
  COUT1(" > [pass #2] reading data");

  vector<int> zoneBcType;

  int current_tri = ss_nst0;
  int nno_count = 0;
  int nfa_count = 0;
  int nozone_count = 0;
  int fazone_count = 0;
  int zoHeader_count = 0;
  int cvzone_count = 0;

  while (fr.advanceToLevel(1)) {

    const int section_header = fr.getNextTokenAsInt(1);

    if (section_header == 0) {}  // comment
    else if (section_header == 1) {}
    else if (section_header == 2) {}
    else if (section_header == 4) {}
    else if (section_header == 37) {}
    else if (section_header == 38) {}
    else if (section_header == 40) {}  // partition info
    else if (section_header == 64) {}
    else if (section_header == 10) {
      // Node Block
      // fomrat: (10 (id start end type [ND]) (x y [z]))
      // [.] indicate optional parameters, depending on dimension
      const int id    = fr.getNextTokenAsHex(2);
      const int start = fr.getNextTokenAsHex(2);
      const int end   = fr.getNextTokenAsHex(2);
      const int type  = fr.getNextTokenAsHex(2);

      if (id != 0) {
        // data present section
        int noND = 3;
        if (fr.getLevel() == 2) noND  = fr.getNextTokenAsHex(2);

        if (noND != 3) {
          WUI(WARN,"Surfer can only import 3D mesh files currently; skipping");
          return -2;
        }

        COUT2("    > node block: id: " << id << ", start: " << start << ", end: " << end << ", type: " << type << ", dimensions: " << noND);

        nno_count += end - start + 1;
        ++nozone_count;

        fr.advanceToLevel(2);

        for (int i = ss_nsp0+start; i <= ss_nsp0+end; ++i) {
          FOR_J3 xsp[i-1][j] = fr.getNextTokenAsDouble(2);
        }
      }
    }
    else if (section_header == 12) {
      // Cell Block
      // format: (12 (id start end type elem-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int type   = fr.getNextTokenAsHex(2);

      if (id == 0) {
        // declaration section, elem-type is omitted
        if (end != 0) {
          COUT2("    > global cell declaration: " << end << " cells");
        }
        else {
          COUT2("    > global cell declaration: boundary mesh (no cell data)");
        }
      }
      else {
        // data present section
        const int elType   = fr.getNextTokenAsHex(2);

        COUT2("    > cell data block: id: " << id << ", start: " << start << ", end: " << end << ", cell-type: " << type << ", element-type: " << elType);

        ++cvzone_count;
      }

    }
    else if (section_header == 13) {
      // Face Block Section
      // format: (13 (id start end bc-type fa-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int type   = fr.getNextTokenAsHex(2);  // bc-type for non-declaration blocks

      if (id == 0) {
        // declaration section, bc-type is omitted
        // set the face count for checking
        assert(start == 1);
        COUT2("    > face declaration: " << end << " faces");
        if (tri_internal) COUT2("    > internal faces skipped: " << tri_internal << " faces");
        COUT2("    > decimated into: " << tri_count << " surface triangles");

        // allocate the surface information
        nst += tri_count;
        growNstData(nst,ss_nst0);
        for (int ist = ss_nst0; ist < nst; ++ist) znost[ist] = -1;

      }
      else {
        // data present section
        const int faType   = fr.getNextTokenAsHex(2);

        COUT2("    > ascii face data block: id: " << id << ", start: " << start << ", end: " << end << ", bc-type: " << type << ", face-type: " << faType);

        bool skipZone = false;
        set<int>::iterator it;
        it = skipFaceZoneIds.find(id);
        if (it!=skipFaceZoneIds.end()) {
          // this zone should be skipped b/c it is internal
          skipZone = true;
        }

        fr.advanceToLevel(2);

        // create unique Zone identifier
        const int izone = zoneVec.size();
        if (!skipZone) {
          nfa_count += end - start + 1;
          ++fazone_count;
          ostringstream oss;
          oss << "X" << id;
          zoneVec.push_back(SurfaceZone(oss.str()));
        }

        if (faType >= 2 && faType <=4) {
          // list all nodes, then two cells
          int x[faType+2];
          for (int ifa = start; ifa <= end; ++ifa) {
            for (int i = 0; i < faType; ++i) {
              x[i] = fr.getNextTokenAsHex(2);
              assert((x[i] >= 1) && (x[i] <= nsp));
            }
            for (int i = faType; i < (faType+2); ++i) {
              x[i] = fr.getNextTokenAsHex(2);
            }

            if (!skipZone) {
              // create the tris for this polygon
              for (int ist = 0; ist < faType-2; ++ist) {
                spost[current_tri][0] = ss_nsp0+x[0]-1;     // faces and nodes are 1-indexed
                spost[current_tri][1] = ss_nsp0+x[ist+1]-1; // faces and nodes are 1-indexed
                spost[current_tri][2] = ss_nsp0+x[ist+2]-1; // faces and nodes are 1-indexed

                znost[current_tri] = izone;
                ++current_tri;
              }
            }
          }
        }
        else if (faType == 0 || faType == 5) {
          // list all nodes, then two cells

          for (int ifa = start; ifa <= end; ++ifa) {
            int no_per_fa = fr.getNextTokenAsHex(2);

            int x[no_per_fa+2];
            for (int i = 0; i < no_per_fa; ++i) {
              x[i] = fr.getNextTokenAsHex(2);
              assert((x[i] >= 1) && (x[i] <= nsp));
            }
            for (int i = no_per_fa; i < (no_per_fa+2); ++i) {
              x[i] = fr.getNextTokenAsHex(2);
            }

            if (!skipZone) {
              // create the tris for this polygon
              for (int ist = 0; ist < no_per_fa-2; ++ist) {
                spost[current_tri][0] = ss_nsp0+x[0]-1;     // faces and nodes are 1-indexed
                spost[current_tri][1] = ss_nsp0+x[ist+1]-1; // faces and nodes are 1-indexed
                spost[current_tri][2] = ss_nsp0+x[ist+2]-1; // faces and nodes are 1-indexed

                znost[current_tri] = izone;
                ++current_tri;
              }
            }
          }
        }
        else {
          WUI(WARN,"unsupported boundary face type detected (faType = " << faType << ")\n");
        }
      }
    }
    else if ((section_header == 2010) || (section_header == 3010)) {
      // binary node block...
      const int id    = fr.getNextTokenAsHex(2);
      const int start = fr.getNextTokenAsHex(2);
      const int end   = fr.getNextTokenAsHex(2);
      const int type  = fr.getNextTokenAsHex(2);
      int noND = 3;
      if (fr.getLevel() == 2) noND  = fr.getNextTokenAsHex(2);

      nno_count += end - start + 1;
      ++nozone_count;

      COUT2("    > binary node data block: id: " << id << ", start: " << start << ", end: " << end << ", type: " << type << ", dimensions: " << noND);

      fr.advanceToLevel(2);

      // ------------------------------------------
      // set surface node stuff...
      if (noND != 3) {
        WUI(WARN,"Surfer can only import 3D mesh files currently; skipping");
        return -2;
      }
      assert((start >= 1) && (end <= (nsp-ss_nsp0))); // ensures there is one node zone -- support multiple in the future
      //assert(nozone_count == 1);
      assert((type == 1) || (type == 2)); // regular or boundary nodes

      if (section_header == 3010) {
        for (int i = ss_nsp0+start; i <= ss_nsp0+end; ++i) fr.readBinaryNodes<double>(xsp[i-1],3); // nodes have 0 indexing, Fluent uses 1-indexing
      }
      else {
        assert(section_header == 2010);
        double x[3];
        for (int i = ss_nsp0+start; i <= ss_nsp0+end; ++i) {
          fr.readBinaryNodes<float>(x,3); // nodes have 0 indexing, Fluent uses 1-indexing
          FOR_J3 xsp[i-1][j] = x[j];
        }
      }
      // ------------------------------------------

    }
    else if ((section_header == 2012) || (section_header == 3012)) {
      // Need to properly stride across binary sections
      // (2012 (id start end bc-type fa-type) ())
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int type   = fr.getNextTokenAsHex(2);
      const int elemType = fr.getNextTokenAsHex(2);
      UNUSED(id);
      UNUSED(type);

      // only need to read cell description if of mixed-type
      if (elemType == 0) {
        fr.advanceToLevel(2);

        int cellType;
        for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(cellType);
        }
      }
      // otherwise cell block is empty

      // ------------------------------------------

    }
    else if ((section_header == 2013) || (section_header == 3013)) {
      // binary faces block...
      // (2013 (id start last bc-type face-type) no n1 n2 c0 c1)
      const int id     = fr.getNextTokenAsHex(2);
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int bcType = fr.getNextTokenAsHex(2);
      const int faType = fr.getNextTokenAsHex(2);

      COUT2("    > binary face data block: id: " << id << ", start: " << start << ", end: " << end << ", bc-type: " << bcType << ", face-type: " << faType);

      fr.advanceToLevel(2);

      bool skipZone = false;
      set<int>::iterator it;
      it = skipFaceZoneIds.find(id);
      if (it!=skipFaceZoneIds.end()) {
        // this zone should be skipped b/c it is internal
        skipZone = true;
      }

      // ----------------------------------------------
      // update surface
      const int izone = zoneVec.size();
      if (!skipZone) {
        nfa_count += end - start + 1;
        ++fazone_count;
        ostringstream oss;
        oss << "X" << id;
        zoneVec.push_back(SurfaceZone(oss.str()));
      }
      // ----------------------------------------------

      if (faType >= 2 && faType <= 4) {
        int x[faType+2];
        for (int ifa = start; ifa <= end; ++ifa) {
          fr.readBinaryFaces(x,faType);

          for (int i = 0; i < faType; ++i) {
            assert((x[i] >= 1) && (x[i] <= nsp));
          }

          if (!skipZone) {
            // create the tris for this polygon
            for (int ist = 0; ist < faType-2; ++ist) {
              spost[current_tri][0] = ss_nsp0+x[0]-1;     // faces and nodes are 1-indexed
              spost[current_tri][1] = ss_nsp0+x[ist+1]-1; // faces and nodes are 1-indexed
              spost[current_tri][2] = ss_nsp0+x[ist+2]-1; // faces and nodes are 1-indexed

              assert(znost[current_tri] == -1);
              znost[current_tri] = izone;
              ++current_tri;
            }
          }
        }
      }
      else if (faType == 0 || faType == 5) {
        // first variable declares number of points per face
        // still only support tris
        int no_per_fa;
        for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(no_per_fa);

          int x[no_per_fa+2];
          fr.readBinaryFaces(x,no_per_fa);
          for (int i = 0; i < no_per_fa; ++i) {
            assert((x[i] >= 1) && (x[i] <= nsp));
          }

          if (!skipZone) {
            // create the tris for this polygon
            for (int ist = 0; ist < no_per_fa-2; ++ist) {
              spost[current_tri][0] = ss_nsp0+x[0]-1;     // faces and nodes are 1-indexed
              spost[current_tri][1] = ss_nsp0+x[ist+1]-1; // faces and nodes are 1-indexed
              spost[current_tri][2] = ss_nsp0+x[ist+2]-1; // faces and nodes are 1-indexed

              assert(znost[current_tri] == -1);
              znost[current_tri] = izone;
              ++current_tri;
            }
          }
        }
      }

    }
    else if (section_header == 2041) {
      // Need to properly stride across binary sections
      // (2041 (start end) ())
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if ((section_header == 2058) || (section_header == 2059)) {
      // Need to properly stride across binary sections
      // (2059/58 (start end parent-zone child-zone) ())
      // where zones are either face or cell blocks, depending
      const int start = fr.getNextTokenAsHex(2);
      const int end = fr.getNextTokenAsHex(2);
      const int parentZone = fr.getNextTokenAsHex(2);
      const int childZone = fr.getNextTokenAsHex(2);
      UNUSED(parentZone);
      UNUSED(childZone);

      fr.advanceToLevel(2);

      int n_children;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(n_children);
        for (int j=0; j<n_children; ++j) fr.getNextBinary<int,int>();  // don't need to store this info
      }
    }
    else if (section_header == 2040) {
      // Need to properly stride across binary sections
      // (2040 (zone-id first-index last-index partition-count) ())
      const int zone_id  = fr.getNextTokenAsHex(2);
      const int start    = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);
      const int part_count    = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if (section_header == 2041) {
      // Need to properly stride across binary sections
      // (2041 (start end) ())
      const int start  = fr.getNextTokenAsHex(2);
      const int end    = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if (section_header == 2018 || section_header == 3018) {
      // Need to properly stride across binary sections
      // (2041 (start end) ())
      const int start = fr.getNextTokenAsHex(2);
      const int end = fr.getNextTokenAsHex(2);
      const int zone_id  = fr.getNextTokenAsHex(2);
      const int shadow_zone_id = fr.getNextTokenAsHex(2);

      fr.advanceToLevel(2);

      // assumes an int flag
      int tmp_flag;
      for (int i = start; i <= end; ++i) {
        // face correspondence between shadow zones
        fr.readBinaryInt(tmp_flag);
        fr.readBinaryInt(tmp_flag);
      }
    }
    else if ((section_header == 45) || (section_header == 39))  {
      // format:
      // (45 (id zone-type zone-name) ())
      const int id     = fr.getNextTokenAsInt(2);
      const string znType = fr.getNextTokenAsString(2);
      const string znName = fr.getNextTokenAsString(2);

      COUT2("    > zone header: zone-id: " << id << ", type: " << znType << ", name: " << znName);

      ostringstream oss;
      oss << "X" << id;
      const string oss_str = oss.str();

      // look for this zone in this->zoneVec...
      int i;
      for (i = 0; i < int(this->zoneVec.size()); ++i) {
        if (this->zoneVec[i].getName() == oss_str) {
          this->zoneVec[i].setName(znName);
          break;
        }
      }
      // make sure we found it...
      //assert(i < this->zoneVec.size()); # assumes sequential ordering from 1...

      if (i< int(this->zoneVec.size())) ++zoHeader_count;

    }
    else {
      WUI(WARN," > skipping unrecognized block id \"" << section_header << "\"\n");
    }
  }

  // Count Checking
  COUT1(" > initFromMsh Summary:");
  COUT1("     " << nno_count << "/" << (nsp-ss_nsp0) << " nodes read");
  COUT1("     " << current_tri << "/" << (nst-ss_nst0) << " surface tris read");
  COUT1("     " << nozone_count << " node zones read");
  COUT1("     " << fazone_count << " face zones read");
  COUT1("     " << cvzone_count << " cell zones read");
  COUT1("     " << zoHeader_count << " zone headers read");

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  // flip all added tris because Fluent convention is opposite of ours
  st_flag.resize(nst);
  st_flag.setAll(0);
  for (int ist = ss_nst0; ist < nst; ++ist) st_flag[ist] = 1;
  flipFlaggedTris();

  return 0;  // success
}
