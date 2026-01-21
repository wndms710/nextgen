#ifndef CAS_DAT_READER_HPP
#define CAS_DAT_READER_HPP
#include "parsers/fluent/CasReader.hpp"

#define MAX_NO_PER_FA      32

class CasDatReader {

private:

  // storage for striped partition
  int8 nno_global;
  int8 nfa_global;
  int8 * noora;
  int8 * faora;
  int nfa;
  int nfa_dup;
  int * noofa_i;
  int8 * noofa_v_global;
  int8 (*noofa_v_global_tmp)[MAX_NO_PER_FA];
  int noofa_s;
  int8 (*cvofa_global)[2];
  int8 * fa_parent_list;
  int nno;
  double (*x_no)[3];

  int force_byte_swap;  // 0: use parsed args, 1: don't swap, 2: swap

  CasReader * crdr;
  int current_data_id; // store current variable data id
  int cv_zones_set;  // keep track of how many cv_zones have been populated
  set<int> cell_block_ids;

public:

  CasDatReader() {
    noora = NULL;
    faora = NULL;
    noofa_i = NULL;
    noofa_v_global_tmp = NULL;
    noofa_v_global = NULL;
    cvofa_global = NULL;
    x_no = NULL;
    nno = 0;
    nno_global = 0;
    nfa_global = 0;
    nfa = 0;
    nfa_dup = 0;
    noofa_s = 0;
    fa_parent_list = NULL ;
    crdr = NULL;

    force_byte_swap = 0;

    current_data_id = 0;
    cv_zones_set = 0;
  }

  ~CasDatReader() {
    DELETE(noora);
    DELETE(faora);
    DELETE(noofa_i);
    DELETE(noofa_v_global_tmp);
    DELETE(noofa_v_global);
    DELETE(cvofa_global);
    DELETE(x_no);
    DELETE(fa_parent_list);

    if (crdr != NULL) delete crdr;
  }

  // references so we can init/alloc them within this fxn
  void initFromCas(double (*&x_cv)[3],double*& delta,int8*& cvora,int& ncv, int8& ncv_global,const string filename) {

    if ((mpi_rank == 0)&&(cti_verbose))
      cout << "CasDatReader::initFromCas() " << filename << endl;

    assert(crdr == NULL);
    crdr = new CasReader(filename);

    if (force_byte_swap == 2) crdr->setByteSwap(true);
    else if (force_byte_swap == 1) crdr->setByteSwap(false);

    while (crdr->advanceToLevel(1)) {
      const int section_header = crdr->getNextTokenAsInt(1);

      if (section_header == 2) {
        processDimension();
      }
      else if (section_header == 4) {
        processByteSwap();
      }
      else if (section_header == 10) {
        processNodes();
      }
      else if ((section_header == 2010) || (section_header == 3010)) {
        const bool b_double = (section_header == 3010) ? true:false;
        processBinaryNodes(b_double);
      }
      else if (section_header == 12) {
        processCells(x_cv,delta,cvora,ncv,ncv_global);
      }
      else if ((section_header == 2012) || (section_header == 3012)) {
        processBinaryCells();
      }
      else if (section_header == 13) {
        processFaces();
      }
      else if ((section_header == 2013) || (section_header == 3013)) {
        processBinaryFaces();
      }
      else if ((section_header == 2018) || (section_header == 3018)) {
        processBinaryShadowZones();
      }
      else if (section_header == 2040) {
        processBinaryCellPartition();
      }
      else if (section_header == 2041) {
        processBinaryFortyOne();
      }
      else if (section_header == 59) {
        processFaceTree();
      }
      else if (section_header == 2059) {
        processBinaryFaceTree();
      }
      else if (section_header == 58) {
        // ignore cell tree
      }
      else if (section_header == 2058) {
        // need to stride across this section
        processBinaryCellTree();
      }
      else if (section_header == 2061) {
        // need to stride across this section
        processBinaryInterfaceParents();
      }
      else {
        COUT1(" > unprocessed Fluent data block with index: " << section_header);
      }
    }

    finalizeAndCheck();

    computeXcvDelta(x_cv,delta,cvora,ncv,ncv_global);

    delete crdr; crdr = NULL;

    if (mpi_rank == 0) cout << " > ncv_global: " << ncv_global << endl;
    MiscUtils::dumpRange(x_cv,ncv,"x_cv from cas/msh");
    MiscUtils::dumpRange(delta,ncv,"delta from cas/msh");
  }

  void initFromDat(list<DnData>& dnList,list<Dn3Data>& dn3List,const string filename,const set<string>& varNameSet,const int8* cvora,const int ncv,const int8 ncv_global) {

    if ((mpi_rank == 0)&&(cti_verbose))
      cout << "CasDatReader::initFromDat() " << filename << endl;

    assert(crdr == NULL);
    crdr = new CasReader(filename);

    while (crdr->advanceToLevel(1)) {
      const int section_header = crdr->getNextTokenAsInt(1);

      if (section_header == 4) {
        processByteSwap();
      }
      else if (section_header == 33) {
        processGridSize(ncv_global);
      }
      else if (section_header == 300) {
        processData(dnList,dn3List,varNameSet,cvora,ncv,ncv_global);
      }
      else if ((section_header == 2300) || (section_header == 3300)) {
        bool read_double = (section_header == 3300) ? true:false;
        processBinaryData(dnList,dn3List,varNameSet,cvora,ncv,ncv_global,read_double);
      }
      else if ((section_header == 2301) || (section_header == 2302) || (section_header == 3301) || (section_header == 3302)) {
        const bool read_double = ((section_header == 3301) || (section_header == 3302)) ? true:false;
        processBinaryResiduals(read_double);
      }
      else {
        COUT1("unrecognized dat-file section with id: " << section_header);
      }
    }

    delete crdr; crdr = NULL;

    MPI_Barrier(mpi_comm);
    COUT1("CasDatReader::initFromDat() done");
  }

  void setForceByteSwap(const int b_swap) {
    force_byte_swap = b_swap;
  }

private:

  void processDimension() {
    const int dim = crdr->getNextTokenAsInt(1);
    COUT2(" > dimensionality: " << dim);

    if (dim != 3) {
      CERR("only 3D meshes are currently supported");
    }

  }

  void processByteSwap() {
    //COUT1("CasDatReader::processByteSwap");

    const int endianFlag = crdr->getNextTokenAsInt(2);
    if (endianFlag != 60) {
      if (force_byte_swap == 0) {
        crdr->setByteSwap(true);
        COUT2(" > byte-swapping binary records");
      }
    }
  }

  void processGridSize(const int8 ncv_global) {
    //COUT1("CasDatReader::processGridSize");

    const int n_cells = crdr->getNextTokenAsInt(2);

    // cell count...
    if (n_cells != ncv_global) {
      CERR("mismatch b/w .cas and .dat cell counts");
    }

    // these later checks are relevant only if we init'ed from .cas/.msh...

    const int n_faces = crdr->getNextTokenAsInt(2);
    // face count...
    if ((nfa_global > 0) && (n_faces != nfa_global)) {
      CERR("mismatch b/w .cas and .dat face counts.");
    }

    const int n_nodes = crdr->getNextTokenAsInt(2);
    // node count...
    if ((nno_global > 0) && (n_nodes != nno_global)) {
      CERR("mismatch b/w .cas and .dat node counts.");
    }

    COUT2(" > grid size: " << n_cells << " (cells), " << n_faces << " (faces), " << n_nodes << " (nodes)");

  }

  // ASCII node section
  void processNodes() {
    //COUT1("CasDatReader::processNodes");

    const int id    = crdr->getNextTokenAsHex(2);
    const int8 start = crdr->getNextTokenAsHex(2);
    const int8 end   = crdr->getNextTokenAsHex(2);
    const int type  = crdr->getNextTokenAsHex(2);

    if (id == 0) {
      assert(start == 1);

      // the global node count is end...
      // we will distribute the reading of the nodes across the processors...

      // just lump all nodes together for now...
      nno_global = end;
      COUT1(" > nno_global: "<< nno_global);

      // nodes are striped as evenly as possible across all processors...
      assert(noora == NULL);
      MiscUtils::calcUniformDist(noora, nno_global, mpi_size);

      assert(nno == 0);
      nno = noora[mpi_rank+1] - noora[mpi_rank];

      // need these to (eventually compute x_cv)...
      assert(x_no == NULL);
      x_no = new double[nno][3];
    }
    else {
      // this block contains data
      assert(id != 0);
      COUT2(" > ascii node block with id: " << id << " (start: " << start << ", end: " << end << ", type: " << type << ")");

      // optionally specified dimensionality
      int noND = 3;
      if (crdr->getLevel() == 2) noND = crdr->getNextTokenAsHex(2);

      if (noND != 3) {
        CERR("2D nodal data is not supported");
      }

      crdr->advanceToLevel(2);

      // next are the global indices of this node range...
      // fluent indices are 1-indexed
      int8 istart = start - 1;
      int8 iend = end - 1;

      // only parse the remainder if we own nodes in this range...
      //if ((istart < noora[mpi_rank+1])&&(iend >= noora[mpi_rank])) {

        assert(x_no != NULL );  // should have been allocated in ID=0 block

        for (int8 ino_global = istart; ino_global <= iend; ino_global++) {
          const int ino = ino_global - noora[mpi_rank];
          if ((ino >= 0)&&(ino < nno)) {
            // parse 3 doubles...
            FOR_I3 x_no[ino][i] = crdr->getNextTokenAsDouble(2);
          }
          else {
            // advance 3 doubles...
            FOR_I3 crdr->getNextTokenAsDouble(2);
          }
        }
      //}
    }
  }

  void processBinaryNodes(const bool b_double) {

    const int id    = crdr->getNextTokenAsHex(2);
    const int8 start = crdr->getNextTokenAsHex(2);
    const int8 end   = crdr->getNextTokenAsHex(2);
    const int type  = crdr->getNextTokenAsHex(2);

    // this block contains data
    assert(id != 0);
    COUT2(" > binary node block with id: " << id << " (start: " << start << ", end: " << end << ", type: " << type << ")");

    // optionally specified dimensionality
    int noND = 3;
    if (crdr->getLevel() == 2) noND = crdr->getNextTokenAsHex(2);

    if (noND != 3) {
      CERR("2D nodal data is not supported");
    }

    crdr->advanceToLevel(2);

    // next are the global indices of this node range...
    // fluent indices are 1-indexed
    int8 istart = start - 1;
    int8 iend = end - 1;

    // only parse the remainder if we own nodes in this range...
    // if ((istart < noora[mpi_rank+1])&&(iend >= noora[mpi_rank])) {
    for (int8 ino_global = istart; ino_global <= iend; ino_global++) {

      assert(x_no != NULL );  // should have been allocated in ID=0 block

      // for (int8 ino_global = istart; ino_global <= iend; ino_global++) {
        const int ino = ino_global - noora[mpi_rank];
        if ((ino >= 0)&&(ino < nno)) {
          // parse 3 doubles...
          if (b_double) crdr->readBinaryNodes<double>(x_no[ino],3);
          else crdr->readBinaryNodes<float>(x_no[ino],3);
        }
        else {
          // advance 3 doubles...
          double tmp[3];
          if (b_double) crdr->readBinaryNodes<double>(tmp,3);
          else crdr->readBinaryNodes<float>(tmp,3);
        }
      // }
    }
  }

  // referenced so we can set/alloc them within the fxn
  void processCells(double (*&x_cv)[3],double*& delta,int8*& cvora,int& ncv,int8& ncv_global) {
    //COUT1("CasDatReader::processCvs");

    const int id     = crdr->getNextTokenAsHex(2);
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int type   = crdr->getNextTokenAsHex(2);

    if (id == 0) {
      // declaration block

      // cell count...
      assert(start == 1);

      ncv_global = end;
      COUT1(" > ncv_global: "<< ncv_global);

      MiscUtils::calcUniformDist(cvora, ncv_global, mpi_size);

      assert(ncv == 0);
      ncv = cvora[mpi_rank+1] - cvora[mpi_rank];

      // allocate x_cv and delta now that we know how many centroids we need
      x_cv = new double[ncv][3];
      delta = new double[ncv];
    }
    else {
      // data block
      assert(id != 0); // never zero
      cell_block_ids.insert(id);
      COUT2(" > ascii cell block with id: " << id << " (start: " << start << ", end: " << end << ", type: " << type << ")");
      // do nothing because we will build cell information strictly from the faces & nodes
    }
  }

  void processBinaryCells() {
    //COUT1("CasDatReader::processCvs");

    // Need to properly stride across binary sections
    // (2012 (id start end bc-type fa-type) ())
    const int id     = crdr->getNextTokenAsHex(2);
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int type   = crdr->getNextTokenAsHex(2);
    const int elemType = crdr->getNextTokenAsHex(2);
    UNUSED(id);
    UNUSED(type);

    cell_block_ids.insert(id);
    COUT2(" > binary cell block with id: " << id << " (start: " << start << ", end: " << end << ", bc-type: " << type << ", element-type: " << elemType << ")");
    // only need to read cell description if of mixed-type
    if (elemType == 0) {
      crdr->advanceToLevel(2);

      int cellType;
      for (int i = start; i <= end; ++i) {
        crdr->readBinaryInt(cellType);
      }
    }
    // otherwise cell block is empty
  }

  void processFaceTree() {
    //COUT1("CasDatReader::processDuplicatedFaces");

    // (59 (start end parentID childID)
    // (nkids kid-id0 kid-id1 ... kid-idn))
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int parentId   = crdr->getNextTokenAsHex(2);
    const int childId   = crdr->getNextTokenAsHex(2);
    UNUSED(parentId);
    UNUSED(childId);

    COUT2(" > ascii face tree section found");

    const int8 istart = start-1;
    const int8 iend   = end-1;

    //TODO this only allows one of these sections to exist...assume will break if multiple due to the reallocation of fa_parent_list...

    if ((istart < faora[mpi_rank+1])&&(iend >= faora[mpi_rank])) {

      int8 is    = max(istart,faora[mpi_rank]);
      int8 ie    = min(iend  ,faora[mpi_rank+1]-1);
      nfa_dup = ie-is+1;
      assert ( nfa_dup > 0 ) ;
      fa_parent_list = new int8[nfa_dup];

      for (int ifa_dup = 0; ifa_dup < nfa_dup; ++ifa_dup)
        fa_parent_list[ifa_dup] = ifa_dup + is ; // global idx ..
    }
    // only use the header to grab redundant faces, ignore body
  }

  void processBinaryFaceTree(const bool b_ignore = false) {
    //COUT1("CasDatReader::processDuplicatedFaces");

    // (2059 (start end parentID childID)
    // (nkids kid-id0 kid-id1 ... kid-idn))
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int parentId   = crdr->getNextTokenAsHex(2);
    const int childId   = crdr->getNextTokenAsHex(2);
    UNUSED(parentId);
    UNUSED(childId);

    const int8 istart = start-1;
    const int8 iend   = end-1;

    //TODO this only allows one of these sections to exist...assume will break if multiple due to the reallocation of fa_parent_list...

    if (!b_ignore) {
      COUT2(" > binary face tree section found");

      if ((istart < faora[mpi_rank+1])&&(iend >= faora[mpi_rank])) {

        int8 is    = max(istart,faora[mpi_rank]);
        int8 ie    = min(iend  ,faora[mpi_rank+1]-1);
        nfa_dup = ie-is+1;
        assert ( nfa_dup > 0 ) ;
        fa_parent_list = new int8[nfa_dup];

        for (int ifa_dup = 0; ifa_dup < nfa_dup; ++ifa_dup)
          fa_parent_list[ifa_dup] = ifa_dup + is ; // global idx ..
      }
    }

    // only use the header to grab redundant faces, stride over body
    crdr->advanceToLevel(2);

    for (int icell = istart; icell <= iend; ++icell) {
      const int n_children = crdr->getNextBinary<int,int>();
      int tmp[n_children];
      crdr->readBinaryInts(&tmp[0],n_children);
    }
  }

  void processBinaryCellTree() {
    COUT2(" > binary cell tree section found");
    processBinaryFaceTree(true);
  }

  void processFaces() {
    //COUT1("CasDatReader::processFaces");

    // format: (13 (id start end bc-type fa-type) ())
    const int id     = crdr->getNextTokenAsHex(2);
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int type   = crdr->getNextTokenAsHex(2);  // bc-type for non-declaration blocks

    if (id == 0) {
      // declaration section
      // face count...
      assert(start == 1);

      // faces are striped as evenly as possible across all processors...
      nfa_global = end;
      COUT1(" > nfa_global: "<< nfa_global);

      assert(faora == NULL);
      MiscUtils::calcUniformDist(faora, nfa_global, mpi_size);

      assert(nfa == 0);
      nfa = faora[mpi_rank+1] - faora[mpi_rank];

      // faced-based connectivity structures...
      assert(noofa_i == NULL); noofa_i = new int[nfa+1];
      assert(cvofa_global == NULL); cvofa_global = new int8[nfa][2];
      assert(noofa_v_global_tmp == NULL); noofa_v_global_tmp = new int8[nfa][MAX_NO_PER_FA];
    }
    else {
      assert(id != 0);

      int8 istart = start - 1;
      int8 iend = end - 1;

      const int faType   = crdr->getNextTokenAsHex(2);
      COUT2(" > ascii face block with id: " << id << " (start: " << start << ", end: " << end << ", bc-type: " << type << ", face-type: " << faType << ")");

      crdr->advanceToLevel(2);

      // only parse the rest if you need to...
      //if ((istart < faora[mpi_rank+1])&&(iend >= faora[mpi_rank])) {

        if ((faType == 0)||(faType == 5)) {
          // this version of the faces has the number
          // of nodes on each line...
          for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

            const int nno_per_face = crdr->getNextTokenAsHex(2);
            if (nno_per_face > MAX_NO_PER_FA) cout << "WARNING! nno_per_face: " << nno_per_face << " exceeds limit of: " << MAX_NO_PER_FA << endl;
            assert(nno_per_face <= MAX_NO_PER_FA);

            int no_list[MAX_NO_PER_FA];
            for (int i = 0; i < nno_per_face; i++) {
              no_list[i] = crdr->getNextTokenAsHex(2) - 1; // use zero-indexing
            }

            // after all faces the cells are listed
            int8 icv0 = crdr->getNextTokenAsHex(2) - 1;
            int8 icv1 = crdr->getNextTokenAsHex(2) - 1;

            // get the local face indexing...
            int ifa = ifa_global - faora[mpi_rank];

            // only use the parsed info if we own these faces
            if ((ifa >= 0)&&(ifa < nfa)) {

              // flip the face right away so that any face with a -1 cv
              // has it in the second position, i.e. cvofa_global[ifa][1]...
              if (icv0 >= 0) {
                // put no count in ifa+1 position for CSR struct
                noofa_i[ifa+1] = nno_per_face;
                // reverse the order of the nodes - msh has left-handed rule...
                for (int i = 0; i < nno_per_face; i++)
                  noofa_v_global_tmp[ifa][i] = no_list[nno_per_face-i-1];
                cvofa_global[ifa][0] = icv0;
                cvofa_global[ifa][1] = icv1;
              }
              else {
                assert(icv1 >= 0);
                // put no count in ifa+1 position for CSR struct
                noofa_i[ifa+1] = nno_per_face;
                // use the current face ordering...
                for (int i = 0; i < nno_per_face; i++)
                  noofa_v_global_tmp[ifa][i] = no_list[i];
                cvofa_global[ifa][0] = icv1;
                cvofa_global[ifa][1] = icv0;
              }
            }
          }
        }
        else {
          // these faces have the same number of nodes...
          assert(faType <= MAX_NO_PER_FA);

          // this version of the faces has the number of nodes on each line...
          for (int8 ifa_global = istart; ifa_global <= iend; ifa_global++) {
            int8 no_list[MAX_NO_PER_FA];
            for (int i = 0; i < faType; i++) {
              no_list[i] = crdr->getNextTokenAsHex(2) - 1; // use zero-indexing
            }

            int8 icv0 = crdr->getNextTokenAsHex(2) - 1;
            int8 icv1 = crdr->getNextTokenAsHex(2) - 1;

            // get the local face indexing...
            int ifa = ifa_global - faora[mpi_rank];
            // only use the parsed info if we own these faces
            if ((ifa >= 0)&&(ifa < nfa)) {

              // flip the face right away so that any face with a -1 cv
              // has it in the second position, i.e. cvofa_global[ifa][1]...
              if (icv0 >= 0) {
                // put no count in ifa+1 position for CSR struct
                noofa_i[ifa+1] = faType;
                // reverse the order of the nodes - msh has left-handed rule...
                for (int i = 0; i < faType; i++)
                  noofa_v_global_tmp[ifa][i] = no_list[faType-i-1];
                cvofa_global[ifa][0] = icv0;
                cvofa_global[ifa][1] = icv1;
              }
              else {
                assert(icv1 >= 0);
                // put no count in ifa+1 position for CSR struct
                noofa_i[ifa+1] = faType;
                // use the current face ordering...
                for (int i = 0; i < faType; i++)
                  noofa_v_global_tmp[ifa][i] = no_list[i];
                cvofa_global[ifa][0] = icv1;
                cvofa_global[ifa][1] = icv0;
              }
            }
          }
        }
      }
    //}
  }

  void processBinaryFaces() {
    //COUT1("CasDatReader::processFaces");

    // format: (2013/3013 (id start end bc-type fa-type) ())
    const int id     = crdr->getNextTokenAsHex(2);
    const int8 start  = crdr->getNextTokenAsHex(2);
    const int8 end    = crdr->getNextTokenAsHex(2);
    const int type   = crdr->getNextTokenAsHex(2);  // bc-type

    assert(id != 0);

    int8 istart = start - 1;
    int8 iend = end - 1;

    const int faType   = crdr->getNextTokenAsHex(2);
    COUT2(" > binary face block with id: " << id << " (start: " << start << ", end: " << end << ", bc-type: " << type << ", face-type: " << faType << ")");

    crdr->advanceToLevel(2);

    // only parse the rest if you need to...
    // if ((istart < faora[mpi_rank+1])&&(iend >= faora[mpi_rank])) {
    for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

      if ((faType == 0)||(faType == 5)) {
        // this version of the faces has the number
        // of nodes on each line...
        // for (int ifa_global = istart; ifa_global <= iend; ifa_global++) {

          const int nno_per_face = crdr->getNextBinary<int,int>();
            if (nno_per_face > MAX_NO_PER_FA) cout << "WARNING! nno_per_face: " << nno_per_face << " exceeds limit of: " << MAX_NO_PER_FA << endl;
          assert(nno_per_face <= MAX_NO_PER_FA);

          int no_list[MAX_NO_PER_FA];
          for (int i = 0; i < nno_per_face; i++) {
            no_list[i] = crdr->getNextBinary<int,int>() - 1; // use zero-indexing
          }

          // after all faces the cells are listed
          int8 icv0 = crdr->getNextBinary<int,int>() - 1;
          int8 icv1 = crdr->getNextBinary<int,int>() - 1;

          // get the local face indexing and update if we own
          int ifa = ifa_global - faora[mpi_rank];
          if ((ifa >= 0)&&(ifa < nfa)) {

            // flip the face right away so that any face with a -1 cv
            // has it in the second position, i.e. cvofa_global[ifa][1]...
            if (icv0 >= 0) {
              // put no count in ifa+1 position for CSR struct
              noofa_i[ifa+1] = nno_per_face;
              // reverse the order of the nodes - msh has left-handed rule...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_global_tmp[ifa][i] = no_list[nno_per_face-i-1];
              cvofa_global[ifa][0] = icv0;
              cvofa_global[ifa][1] = icv1;
            }
            else {
              assert(icv1 >= 0);
              // put no count in ifa+1 position for CSR struct
              noofa_i[ifa+1] = nno_per_face;
              // use the current face ordering...
              for (int i = 0; i < nno_per_face; i++)
                noofa_v_global_tmp[ifa][i] = no_list[i];
              cvofa_global[ifa][0] = icv1;
              cvofa_global[ifa][1] = icv0;
            }
          }
        // }
      }
      else {
        // these faces have the same number of nodes...
        assert(faType <= MAX_NO_PER_FA);

        // this version of the faces has the number of nodes on each line...
        // for (int8 ifa_global = istart; ifa_global <= iend; ifa_global++) {
          int8 no_list[MAX_NO_PER_FA];
          for (int i = 0; i < faType; i++) {
            no_list[i] = crdr->getNextBinary<int,int>() - 1; // use zero-indexing
          }

          int8 icv0 = crdr->getNextBinary<int,int>() - 1;
          int8 icv1 = crdr->getNextBinary<int,int>() - 1;

          // get the local face indexing and update if we own
          int ifa = ifa_global - faora[mpi_rank];
          if ((ifa >= 0)&&(ifa < nfa)) {

            // flip the face right away so that any face with a -1 cv
            // has it in the second position, i.e. cvofa_global[ifa][1]...
            if (icv0 >= 0) {
              // put no count in ifa+1 position for CSR struct
              noofa_i[ifa+1] = faType;
              // reverse the order of the nodes - msh has left-handed rule...
              for (int i = 0; i < faType; i++)
                noofa_v_global_tmp[ifa][i] = no_list[faType-i-1];
              cvofa_global[ifa][0] = icv0;
              cvofa_global[ifa][1] = icv1;
            }
            else {
              assert(icv1 >= 0);
              // put no count in ifa+1 position for CSR struct
              noofa_i[ifa+1] = faType;
              // use the current face ordering...
              for (int i = 0; i < faType; i++)
                noofa_v_global_tmp[ifa][i] = no_list[i];
              cvofa_global[ifa][0] = icv1;
              cvofa_global[ifa][1] = icv0;
            }
          }
        }
      // }
    }
  }

  void processBinaryCellPartition() {
    // (2040 (zone-id first-index last-index partition-count) ())
    const int zone_id  = crdr->getNextTokenAsHex(2);
    const int start    = crdr->getNextTokenAsHex(2);
    const int end    = crdr->getNextTokenAsHex(2);
    const int part_count    = crdr->getNextTokenAsHex(2);


    COUT2(" > binary cell partition section found");

    crdr->advanceToLevel(2);

    // assumes an int flag
    int tmp_flag;
    for (int i = start; i <= end; ++i) {
      crdr->readBinaryInt(tmp_flag);
    }
  }

  void processBinaryShadowZones() {
    // Need to properly stride across binary sections
    const int start = crdr->getNextTokenAsHex(2);
    const int end = crdr->getNextTokenAsHex(2);
    const int zone_id  = crdr->getNextTokenAsHex(2);
    const int shadow_zone_id = crdr->getNextTokenAsHex(2);

    crdr->advanceToLevel(2);
    COUT2(" > binary \"18\" section found");

    // assumes an int flag
    int tmp_flag;
    for (int i = start; i <= end; ++i) {
      // face correspondence between shadow zones
      crdr->readBinaryInt(tmp_flag);
      crdr->readBinaryInt(tmp_flag);
    }
  }

  void processBinaryFortyOne() {
    // (2041 (start end) ())
    const int start  = crdr->getNextTokenAsHex(2);
    const int end    = crdr->getNextTokenAsHex(2);

    COUT2(" > binary \"41\" section found");

    crdr->advanceToLevel(2);

    // assumes an int flag
    int tmp_flag;
    for (int i = start; i <= end; ++i) {
      crdr->readBinaryInt(tmp_flag);
    }
  }

  void processBinaryInterfaceParents() {
    // (2061 (face-id0 face-id1) (
    //  parent-id0 parent-id1
    // ))
    const int start  = crdr->getNextTokenAsHex(2);
    const int end    = crdr->getNextTokenAsHex(2);

    COUT2(" > binary interface parents block found");

    crdr->advanceToLevel(2);

    // assumes an int flag
    int tmp_flag[2];
    for (int i = start; i <= end; ++i) {
      crdr->readBinaryInts(tmp_flag,2);
    }
  }

  void processData(list<DnData>& dnList,list<Dn3Data>& dn3List,const set<string>& varNameSet,const int8* cvora,const int ncv,const int8 ncv_global) {
    COUT1("processData()");

    int split3_index = -1; // index for vector data spread over 3 scalars (u[3] -> u,v,w)
    bool found = false; // are we good to start reading

    const int variable_id = crdr->getNextTokenAsInt(2);
    const int zone_id = crdr->getNextTokenAsInt(2);
    const int var_size = crdr->getNextTokenAsInt(2);
    const int time_level = crdr->getNextTokenAsInt(2);
    const int n_phases = crdr->getNextTokenAsInt(2);
    const int8 start = crdr->getNextTokenAsInt(2);
    const int8 end = crdr->getNextTokenAsInt(2);

    COUT2(" > ascii data found for variable-id: " << variable_id << " (zone-id: " << zone_id << ", var-size: " << var_size << ", start: " << start << ", end: " << end << ")");

    crdr->advanceToLevel(2);

    if (cell_block_ids.find(zone_id) != cell_block_ids.end()) {
      if (current_data_id != variable_id) {
        // reset all cell bock flags
        current_data_id = variable_id;
        cv_zones_set = 0;
      }
      else {
        ++cv_zones_set;  // assumes there won't be any redundant blocks...
      }

      // pressure
      if (variable_id == 1) {
        // requested/registered?
        if (varNameSet.find("p") != varNameSet.end()) {
          if (dnList.back().name == "p") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("p"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
      // velocity
      else if (variable_id == 2) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          if (dn3List.back().name == "u") {
            // assume for now that variable cell blocks are contiguous...
            assert(dn3List.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dn3List.push_back(Dn3Data("u"));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[ncv][3];
          }
          found = true;
        }
      }
      // x-velocity
      else if (variable_id == 111) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          // we are assuming that x-velocity is listed first (all cell blocks), then y-velocity, then z-velocity
          if (dn3List.back().name == "u") {
            // assume for now that variable cell blocks are contiguous...
            assert(dn3List.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dn3List.push_back(Dn3Data("u"));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[ncv][3];
          }
          found = true;
          split3_index = 0; // x component
        }
      }
      // y-velocity
      else if (variable_id == 112) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          assert(dn3List.back().name == "u"); // make sure we have u
          found = true;
          split3_index = 1; // y component
        }
      }
      // z-velocity
      else if (variable_id == 113) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          assert(dn3List.back().name == "u"); // make sure we have u
          found = true;
          split3_index = 2; // z component
        }
      }
      // temperature
      else if (variable_id == 3) {
        // requested/registered?
        if (varNameSet.find("T") != varNameSet.end()) {
          if (dnList.back().name == "T") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("T"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
      // density
      else if (variable_id == 101) {
        // requested/registered?
        if (varNameSet.find("rho") != varNameSet.end()) {
          if (dnList.back().name == "rho") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("rho"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
    }
    // mixture fraction
    else if (variable_id == 16) {
      // requested/registered?
      if (varNameSet.find("Z") != varNameSet.end()) {
        if (dnList.back().name == "Z") {
          // assume for now that variable cell blocks are contiguous...
          assert(dnList.back().data != NULL);
        }
        else {
          // not registered yet, so add scalar
          dnList.push_back(DnData("Z"));
          assert(dnList.back().data == NULL);
          dnList.back().data = new double[ncv];
        }
        found = true;
      }
    }
    // progress variable: commented out b/c not available in *.ip reader
    // else if (variable_id == 42) {
    //   // requested/registered?
    //   if (varNameSet.find("C") != varNameSet.end()) {
    //     if (dnList.back().name == "C") {
    //       // assume for now that variable cell blocks are contiguous...
    //       assert(dnList.back().data != NULL);
    //     }
    //     else {
    //       // not registered yet, so add scalar
    //       dnList.push_back(DnData("C"));
    //       assert(dnList.back().data == NULL);
    //       dnList.back().data = new double[ncv];
    //     }
    //     found = true;
    //   }
    // }

    if (found) {
      COUT2("    > storing data; valid cv-zone and variable");

      // size
      if (!((var_size == 1)||(var_size == 3))) {
        cerr << "Error: expect variable size = 1 or 3" << endl;
        throw(-1);
      }

      int8 istart = start - 1;
      int8 iend = end - 1;

      // everybody parses, but only save the data if you own it

      if (var_size == 1) {
        if (split3_index == -1) {
          for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
            const int icv = icv_global - cvora[mpi_rank];
            const double tmp = crdr->getNextTokenAsDouble(2);
            if ((icv >= 0)&&(icv < ncv)) {
              // only store if we own this cv
              dnList.back().data[icv] = tmp;
            }
          }
        }
        else {
          for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
            const int icv = icv_global - cvora[mpi_rank];
            const double tmp = crdr->getNextTokenAsDouble(2);
            if ((icv >= 0)&&(icv < ncv)) {
              // store if we own this cv
              dn3List.back().data[icv][split3_index] = tmp;
            }
          }
        }
      }
      else if (var_size == 3) {
        for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
          const int icv = icv_global - cvora[mpi_rank];
          double tmp[3];
          FOR_I3 tmp[i] = crdr->getNextTokenAsDouble(2);  // everybody reads

          if ((icv >= 0)&&(icv < ncv)) {
            // store if we own this cv
            FOR_I3 dn3List.back().data[icv][i] = tmp[icv];
          }
        }
      }

      MPI_Barrier(mpi_comm);
      // range reporting, but make sure all cv blocks have been processed
      if (cv_zones_set == int(cell_block_ids.size())) {
        if (var_size == 1) {
          if (split3_index == -1) {
            MiscUtils::dumpRange(dnList.back().data,ncv,dnList.back().name);
          }
          else if (split3_index == 2) {
            MiscUtils::dumpRange(dn3List.back().data,ncv,dn3List.back().name);
          }
        }
        else {
          assert(var_size == 3);
          MiscUtils::dumpRange(dn3List.back().data,ncv,dn3List.back().name);
        }
      }
    }
  }

  void processBinaryData(list<DnData>& dnList,list<Dn3Data>& dn3List,const set<string>& varNameSet,const int8* cvora,const int ncv,const int8 ncv_global,const bool b_double) {
    COUT1("CasDatReader::processBinaryData");

    int split3_index = -1; // index for vector data spread over 3 scalars (u[3] -> u,v,w)
    bool found = false; // are we good to start reading

    const int variable_id = crdr->getNextTokenAsInt(2);
    const int zone_id = crdr->getNextTokenAsInt(2);  // could be fa or cv zone
    const int var_size = crdr->getNextTokenAsInt(2);
    const int time_level = crdr->getNextTokenAsInt(2);
    const int n_phases = crdr->getNextTokenAsInt(2);
    const int8 start = crdr->getNextTokenAsInt(2);
    const int8 end = crdr->getNextTokenAsInt(2);

    COUT2(" > binary data found for variable-id: " << variable_id << " (zone-id: " << zone_id << ", var-size: " << var_size << ", start: " << start << ", end: " << end << ")");

    crdr->advanceToLevel(2);

    if (cell_block_ids.find(zone_id) != cell_block_ids.end()) {
      if (current_data_id != variable_id) {
        // reset all cell bock flags
        current_data_id = variable_id;
        cv_zones_set = 0;
      }
      else {
        ++cv_zones_set;  // assumes there won't be any redundant blocks...
      }
      // this data appleis to one of the cv zones we are processing, so read the data for the variables we want to read

      // pressure
      if (variable_id == 1) {
        // requested/registered?
        if (varNameSet.find("p") != varNameSet.end()) {
          if (dnList.back().name == "p") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("p"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
      // velocity
      else if (variable_id == 2) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          if (dn3List.back().name == "u") {
            // assume for now that variable cell blocks are contiguous...
            assert(dn3List.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dn3List.push_back(Dn3Data("u"));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[ncv][3];
          }
          found = true;
        }
      }
      // x-velocity
      else if (variable_id == 111) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          if (dn3List.back().name == "u") {
            // assume for now that variable cell blocks are contiguous...
            assert(dn3List.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dn3List.push_back(Dn3Data("u"));
            assert(dn3List.back().data == NULL);
            dn3List.back().data = new double[ncv][3];
          }
          found = true;
          split3_index = 0; // x component
        }
      }
      // y-velocity
      else if (variable_id == 112) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          assert(dn3List.back().name == "u"); // make sure we have u
          found = true;
          split3_index = 1; // y component
        }
      }
      // z-velocity
      else if (variable_id == 113) {
        // requested/registered?
        if (varNameSet.find("u") != varNameSet.end()) {
          assert(dn3List.back().name == "u"); // make sure we have u
          found = true;
          split3_index = 2; // z component
        }
      }
      // temperature
      else if (variable_id == 3) {
        // requested/registered?
        if (varNameSet.find("T") != varNameSet.end()) {
          if (dnList.back().name == "T") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("T"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
      // density
      else if (variable_id == 101) {
        // requested/registered?
        if (varNameSet.find("rho") != varNameSet.end()) {
          if (dnList.back().name == "rho") {
            // assume for now that variable cell blocks are contiguous...
            assert(dnList.back().data != NULL);
          }
          else {
            // not registered yet, so add scalar
            dnList.push_back(DnData("rho"));
            assert(dnList.back().data == NULL);
            dnList.back().data = new double[ncv];
          }
          found = true;
        }
      }
    }

    if (found) {
      COUT2("    > storing data; valid cv-zone and variable");
    }


    // size
    if (!((var_size == 1)||(var_size == 3))) {
      cerr << "Error: expect variable size = 1 or 3" << endl;
      throw(-1);
    }

    // for binary data we have to stride across it regardless of whether we want to store it or not
    // only save the data to the appropriate buffer when appropriate, but always read

    int8 istart = start - 1;
    int8 iend = end - 1;

    if (var_size == 1) {
      if (split3_index == -1) {
        for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
          const int icv = icv_global - cvora[mpi_rank];
          const double tmp = (b_double) ? crdr->getNextBinary<double,double>():crdr->getNextBinary<double,float>();
          if (found&&(icv >= 0)&&(icv < ncv)) {
            // only store if we own this cv
            dnList.back().data[icv] = tmp;
          }
        }
      }
      else {
        for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
          const int icv = icv_global - cvora[mpi_rank];
          const double tmp = (b_double) ? crdr->getNextBinary<double,double>():crdr->getNextBinary<double,float>();
          if (found&&(icv >= 0)&&(icv < ncv)) {
            // store if we own this cv
            dn3List.back().data[icv][split3_index] = tmp;
          }
        }
      }
    }
    else if (var_size == 3) {
      for (int8 icv_global = istart; icv_global <= iend; icv_global++) {
        const int icv = icv_global - cvora[mpi_rank];
        double tmp[3];
        FOR_I3 tmp[i] = (b_double) ? crdr->getNextBinary<double,double>():crdr->getNextBinary<double,float>();

        if (found&&(icv >= 0)&&(icv < ncv)) {
          // store if we own this cv
          FOR_I3 dn3List.back().data[icv][i] = tmp[icv];
        }
      }
    }

    MPI_Barrier(mpi_comm);
    // range reporting, but make sure all cv blocks have been processed
    if (cv_zones_set == int(cell_block_ids.size())) {
      if (var_size == 1) {
        if (split3_index == -1) {
          MiscUtils::dumpRange(dnList.back().data,ncv,dnList.back().name);
        }
        else if (split3_index == 2) {
          MiscUtils::dumpRange(dn3List.back().data,ncv,dn3List.back().name);
        }
      }
      else {
        assert(var_size == 3);
        MiscUtils::dumpRange(dn3List.back().data,ncv,dn3List.back().name);
      }
    }
  }

  void processBinaryResiduals(const bool b_double) {
    // (2301/2302/3301/3302 (n residual-zone-id size) (
    //  r1
    //  ...
    //  rn
    // ))
    const int n_residuals  = crdr->getNextTokenAsHex(2);
    const int section_id = crdr->getNextTokenAsHex(2);
    const int size = crdr->getNextTokenAsHex(2);

    COUT2(" > binary residuals section found");

    crdr->advanceToLevel(2);

    if (b_double) {
      for (int i = 1; i <= n_residuals; ++i) {
        crdr->getNextBinary<double,double>();
      }
    }
    else {
      for (int i = 1; i <= n_residuals; ++i) {
        crdr->getNextBinary<double,float>();
      }
    }
  }

  void finalizeAndCheck() {
    //COUT1("CasDatReader::finalizeAndCheck");

    // the noofa_i/v CSR structure can now be built...

    noofa_i[0] = 0;
    for (int ifa = 0; ifa < nfa; ifa++) {
      noofa_i[ifa+1] += noofa_i[ifa];
    }
    noofa_s = noofa_i[nfa];
    assert(noofa_v_global == NULL); noofa_v_global = new int8[noofa_s];
    for (int ifa = 0; ifa < nfa; ifa++) {
      int nof_f = noofa_i[ifa];
      int nof_l = noofa_i[ifa+1]-1;
      for (int nof = nof_f; nof <= nof_l; nof++) {
        noofa_v_global[nof] = noofa_v_global_tmp[ifa][nof-nof_f];
      }
    }
    delete[] noofa_v_global_tmp; noofa_v_global_tmp = NULL;

    // dump some information about faces of interest
    if ( nfa_dup > 0 ) {

      int * fa_flag = new int[nfa];

      int noofa_s_tmp = 0;
      int * noofa_i_tmp = NULL ;
      int8 * noofa_v_global_tmp = NULL;
      int8 (*cvofa_global_tmp)[2] = NULL ;

      int nfa_tmp = 0;
      for (int ifa = 0; ifa < nfa ; ++ifa) fa_flag[ifa] = -1;

      for (int ifa_dup = 0 ; ifa_dup < nfa_dup ; ++ifa_dup) {
        const int ifa = fa_parent_list[ifa_dup] - faora[mpi_rank];
        fa_flag[ifa] = ifa;
      }
      delete[] fa_parent_list; fa_parent_list = NULL;
      delete[] faora; faora = NULL;

      for (int ifa = 0 ; ifa < nfa ; ++ifa) {
        if ( fa_flag[ifa] == -1) {
          fa_flag[ifa] = nfa_tmp++;
          noofa_s_tmp += noofa_i[ifa+1]- noofa_i[ifa];
        }
        else {
          fa_flag[ifa] = -1;
        }
      }

      //cout << " New face count = " << nfa_tmp << endl ;

      noofa_i_tmp = new int[nfa_tmp+1];
      noofa_v_global_tmp = new int8[noofa_s_tmp];
      cvofa_global_tmp = new int8[nfa_tmp][2];

      for (int ifa = 0; ifa < nfa_tmp+1; ++ifa)
        noofa_i_tmp[ifa] = -1;

      // now populate everyone...
      noofa_i_tmp[0] = 0;
      for (int ifa = 0; ifa < nfa; ++ifa) {
        if ( fa_flag[ifa] >= 0 ) { // still a valid face
          int ifa_new  = fa_flag[ifa];
          assert(noofa_i_tmp[ifa_new] >= 0);
          noofa_i_tmp[ifa_new+1] = noofa_i_tmp[ifa_new] + noofa_i[ifa+1] - noofa_i[ifa];
          int nof_new  = noofa_i_tmp[ifa_new];
          for (int nof = noofa_i[ifa]; nof < noofa_i[ifa+1]; ++nof) {
            const int j = nof - noofa_i[ifa];
            noofa_v_global_tmp[nof_new+j] = noofa_v_global[nof];
          }
          cvofa_global_tmp[ifa_new][0] = cvofa_global[ifa][0];
          cvofa_global_tmp[ifa_new][1] = cvofa_global[ifa][1];
        }
      }//ifa

      // now reset the face data
      assert( noofa_i_tmp[nfa_tmp] == noofa_s_tmp);
      for (int ifa = 0; ifa < nfa_tmp; ++ifa)
        assert( noofa_i_tmp[ifa+1] >= 0 );

      nfa = nfa_tmp ;
      noofa_s = noofa_s_tmp ;

      delete[] noofa_i; noofa_i = new int[nfa_tmp+1];
      delete[] noofa_v_global; noofa_v_global = new int8[noofa_s_tmp];
      delete[] cvofa_global; cvofa_global = new int8[nfa_tmp][2];

      noofa_i[0] = noofa_i_tmp[0];
      for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
        noofa_i[ifa+1] = noofa_i_tmp[ifa+1];
        cvofa_global[ifa][0] = cvofa_global_tmp[ifa][0];
        cvofa_global[ifa][1] = cvofa_global_tmp[ifa][1];
      }

      for (int i=0 ; i < noofa_s_tmp; ++i)
        noofa_v_global[i] = noofa_v_global_tmp[i];

      delete[] noofa_i_tmp;
      delete[] noofa_v_global_tmp;
      delete[] cvofa_global_tmp;
      delete[] fa_flag;

    }
  }

  void computeXcvDelta(double (*x_cv)[3],double* delta,int8* cvora,int ncv, int8 ncv_global) {
    //COUT1("CasDatReader::computeXcvDelta");

    // ======================================
    // step 1:
    // get the x_no's to fa's striped dist
    // compute area-weight x_fa
    // ======================================

    // get our global nodes
    map<const int8,int> nodeMap;
    int nno_tmp = 0;
    FOR_IFA {
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int8 ino_global = noofa_v_global[nof]; assert((ino_global >= 0)&&(ino_global < nno_global));
        map<const int8, int>::iterator it = nodeMap.find(ino_global);
        if (it == nodeMap.end()) {
          nodeMap[ino_global] = nno_tmp++;
        }
      }
    }
    int8 *ino_global_tmp = new int8[nno_tmp];
    for (map<const int8,int>::const_iterator cit = nodeMap.begin(); cit != nodeMap.end(); ++cit) {
      ino_global_tmp[cit->second] = cit->first;
    }
    int *noofa_v_tmp = new int[noofa_s];
    FOR_IFA {
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int8 ino_global = noofa_v_global[nof]; assert((ino_global >= 0)&&(ino_global < nno_global));
        map<const int8, int>::const_iterator cit = nodeMap.find(ino_global);
        assert(cit != nodeMap.end());
        noofa_v_tmp[nof] = cit->second; // ino_tmp
      }
    }
    nodeMap.clear();

    // use dde to get x_no on the face striping...
    double (*x_no_tmp)[3] = new double[nno_tmp][3];
    {
      DistributedDataExchanger dde(ino_global_tmp,nno_tmp,noora);
      dde.pull(x_no_tmp,x_no);
    }
    delete[] ino_global_tmp;
    delete[] x_no; x_no = NULL;
    delete[] noora; noora = NULL;

    double (*x_fa)[3] = new double[nfa][3];
    double *A_fa = new double[nfa];
    for (int ifa = 0; ifa < nfa; ++ifa) {
      // get face centroid...
      double l_fa = 0.0;
      FOR_I3 x_fa[ifa][i] = 0.0;
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        int nof2 = nof+1;
        if (nof2 == noofa_i[ifa+1]) nof2 = noofa_i[ifa];
        const int ino = noofa_v_tmp[nof]; assert((ino >= 0)&&(ino < nno_tmp));
        const int ino2 = noofa_v_tmp[nof2]; assert((ino2 >= 0)&&(ino2 < nno_tmp));
        const double dl = DIST(x_no_tmp[ino2],x_no_tmp[ino]);
        l_fa += dl;
        FOR_I3 x_fa[ifa][i] += dl*(x_no_tmp[ino][i] + x_no_tmp[ino2][i]);
      }
      FOR_I3 x_fa[ifa][i] /= (2.0*l_fa);
      // get face area
      A_fa[ifa] = 0.0;
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        int nof2 = nof+1;
        if (nof2 == noofa_i[ifa+1]) nof2 = noofa_i[ifa];
        const int ino = noofa_v_tmp[nof]; assert((ino >= 0)&&(ino < nno_tmp));
        const int ino2 = noofa_v_tmp[nof2]; assert((ino2 >= 0)&&(ino2 < nno_tmp));
        const double n_fa[3] = TRI_NORMAL_2(x_fa[ifa],x_no_tmp[ino2],x_no_tmp[ino]);
        A_fa[ifa] += MAG(n_fa);
      }
      A_fa[ifa] *= 0.5;
    }
    delete[] noofa_v_tmp; noofa_v_tmp = NULL;
    delete[] noofa_i; noofa_i = NULL;
    delete[] noofa_v_global; noofa_v_global = NULL;
    delete[] x_no_tmp;

    // ======================================
    // step 2:
    // get the x_fa's to cv's striped dist
    // compute x_cv
    // ======================================

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;

    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int8 icv0_global = cvofa_global[ifa][0];
      int rank0 = -1;
      if (icv0_global >= 0)
        rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
      const int8 icv1_global = cvofa_global[ifa][1];
      int rank1 = -1;
      if (icv1_global >= 0)
        rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
      if (rank0 == rank1) {
	if (rank0 >= 0) send_count[rank0] += 2;
      }
      else {
	if (rank0 >= 0) send_count[rank0] += 2; // sending both anyways to keep bookeeping simple
	if (rank1 >= 0) send_count[rank1] += 2;
      }
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    int8* send_buf_int8 = new int8[send_count_sum];
    assert(send_count_sum%2 == 0);
    double *send_buf_double = new double[send_count_sum*2]; // x_fa[3],A_fa

    FOR_IFA {
      const int8 icv0_global = cvofa_global[ifa][0];
      int rank0 = -1;
      if (icv0_global >= 0)
        rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
      const int8 icv1_global = cvofa_global[ifa][1];
      int rank1 = -1;
      if (icv1_global >= 0)
        rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
      if (rank0 == rank1) {
        if (rank0 >= 0) {
          send_buf_int8[send_disp[rank0]  ] = icv0_global;
          send_buf_int8[send_disp[rank0]+1] = icv1_global;
          FOR_I3 send_buf_double[send_disp[rank0]*2+i] = x_fa[ifa][i];
          send_buf_double[send_disp[rank0]*2+3] = A_fa[ifa];
          send_disp[rank0] += 2;
        }
      }
      else {
        if (rank0 >= 0) {
          send_buf_int8[send_disp[rank0]  ] = icv0_global;
          send_buf_int8[send_disp[rank0]+1] = icv1_global;
          FOR_I3 send_buf_double[send_disp[rank0]*2+i] = x_fa[ifa][i];
          send_buf_double[send_disp[rank0]*2+3] = A_fa[ifa];
          send_disp[rank0] += 2;
        }
        if (rank1 >= 0) {
          send_buf_int8[send_disp[rank1]  ] = icv1_global;
          send_buf_int8[send_disp[rank1]+1] = icv0_global;
          FOR_I3 send_buf_double[send_disp[rank1]*2+i] = x_fa[ifa][i];
          send_buf_double[send_disp[rank1]*2+3] = A_fa[ifa];
          send_disp[rank1] += 2;
        }
      }
    }
    delete[] x_fa;
    delete[] A_fa;

    // rewind

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    // exchange...

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    int8 * recv_buf_int8 = new int8[recv_count_sum];
    MPI_Alltoallv(send_buf_int8,send_count,send_disp,MPI_INT8,
		  recv_buf_int8,recv_count,recv_disp,MPI_INT8,mpi_comm);
    delete[] send_buf_int8;

    FOR_RANK {
      send_count[rank] *= 2;
      send_disp[rank] *= 2;
      recv_count[rank] *= 2;
      recv_disp[rank] *= 2;
    }

    double *recv_buf_double = new double[recv_count_sum*2];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
		  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double;
    delete[] send_disp;
    delete[] send_count;
    delete[] recv_disp;
    delete[] recv_count;

    double *A_cv = new double[ncv];
    FOR_ICV {
      FOR_I3 x_cv[icv][i] = 0.0;
      A_cv[icv] = 0.0;
    }
    const int nfa_tmp = recv_count_sum/2;
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      double x_fa_tmp[3]; FOR_I3 x_fa_tmp[i] = recv_buf_double[ifa*4+i];
      const double A_fa_tmp = recv_buf_double[ifa*4+3];
      const int8 icv0_global = recv_buf_int8[ifa*2];
      if (icv0_global >= 0) {
        const int rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
        if (rank0 == mpi_rank) {
          const int icv0 = icv0_global-cvora[mpi_rank];
          assert(icv0 >= 0 && icv0 < ncv);
          FOR_I3 x_cv[icv0][i] += A_fa_tmp*x_fa_tmp[i];
          A_cv[icv0] += A_fa_tmp;
        }
      }
      const int8 icv1_global = recv_buf_int8[ifa*2+1];
      if (icv1_global >= 0) {
        const int rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
        if (rank1 == mpi_rank) {
          const int icv1 = icv1_global-cvora[mpi_rank];
          assert(icv1 >= 0 && icv1 < ncv);
          FOR_I3 x_cv[icv1][i] += A_fa_tmp*x_fa_tmp[i];
          A_cv[icv1] += A_fa_tmp;
        }
      }
    }

    FOR_ICV {
      assert(A_cv[icv] > 0.0);
      FOR_I3 x_cv[icv][i] /= A_cv[icv];
    }
    delete[] A_cv;

    // go through again to compute delta
    assert(delta != NULL);
    FOR_ICV delta[icv] = 0.0;
    for (int ifa = 0; ifa < nfa_tmp; ++ifa) {
      double x_fa_tmp[3]; FOR_I3 x_fa_tmp[i] = recv_buf_double[ifa*4+i];
      const int8 icv0_global = recv_buf_int8[ifa*2];
      if (icv0_global >= 0) {
        const int rank0 = MiscUtils::getRankInXora(icv0_global,cvora);
        if (rank0 == mpi_rank) {
          const int icv0 = icv0_global-cvora[mpi_rank];
          assert(icv0 >= 0 && icv0 < ncv);
          delta[icv0] = max(delta[icv0],DIST(x_cv[icv0],x_fa_tmp));
        }
      }
      const int8 icv1_global = recv_buf_int8[ifa*2+1];
      if (icv1_global >= 0) {
        const int rank1 = MiscUtils::getRankInXora(icv1_global,cvora);
        if (rank1 == mpi_rank) {
          const int icv1 = icv1_global-cvora[mpi_rank];
          assert(icv1 >= 0 && icv1 < ncv);
          delta[icv1] = max(delta[icv1],DIST(x_cv[icv1],x_fa_tmp));
        }
      }
    }
    delete[] recv_buf_int8;
    delete[] recv_buf_double;
    FOR_ICV {
      assert(delta[icv] > 0.0);
      delta[icv] *= 2.0; // delta is taken as twice the max distance between a face and the cell centroids
    }

    //char dummy[128];
    //sprintf(dummy,"centroids.%02d.dat",mpi_rank);
    //FILE * fp2 = fopen(dummy,"w");
    //FOR_ICV {
    //  fprintf(fp2,"%18.16e %18.16e %18.16e %18.16e\n",x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],delta[icv]);
    //}
    //fclose(fp2);

  }

};

#endif // CAS_DAT_READER_HPP
