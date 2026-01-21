#ifndef _TEOST_HPP_
#define _TEOST_HPP_

#include "CTI.hpp"
#include "DoubleVertex.hpp"

// 64 bits total:

// | bits (2) | edge (2) | value (60) |
#define BITS_OFFSET 62
#define EDGE_OFFSET 60
#define MASK_VALUE 1152921504606846975ll  // 60 bits (2^BITS - 1)
#define MASK_BITS 3  // 2 bits
#define MASK_EDGE 3  // 2 bits
#define MASK_ME_ORIENT 3  // 2 bits


class Teost {
private:
  int nst;  // current size
  uint8 (*teost_data)[3]; // hold the tri/edge index of a matched edge, or -1 (or -ve open edge loop index - 1)
  vector<pair<int,int> > noome_vec; // node-of-multiedge
  vector<int> teome_i_vec; // tri/edge-of-multiedge index
  vector<uint8> teome_v_vec; // tri/edge-of-multiedge value

  // typedef pair<DoubleVertex,DoubleVertex> edgeNodes;
  // set<edgeNodes> nooie_set; // nodes-of-imprinted-edge

public:
  bool b_non_manifold[3];  // non-manifold topology from edge data

  Teost() {
    nst = 0;
    teost_data = NULL;
    FOR_I3 b_non_manifold[i] = false;
  }

  ~Teost() {
    clear();
  }

  void clear() {
    nst = 0;
    DELETE(teost_data);
    FOR_I3 b_non_manifold[i] = false;

    noome_vec.clear();
    teome_i_vec.clear();
    teome_v_vec.clear();
  }

  // bits == 0: open edge
  // bits == 1: normal face connection (i.e. tri normals are properly aligned)
  // bits == 2: flipped face connection (i.e. tri normals are misaligned)
  // bits == 3: index into multiedge structure

  bool build(const int nsp,const int _nst,const int (*spost)[3]) {
    assert(teost_data == NULL);
    assert(noome_vec.empty());
    assert(teome_i_vec.empty());
    assert(teome_v_vec.empty());

    nst = _nst;
    teost_data = new uint8[nst][3];
    teome_i_vec.resize(1); // tri/edge-of-multiedge index
    teome_i_vec[0] = 0;

    int * neosp_i = new int[nsp+1]; // node,edge-of-surface-point
    for (int isp = 0; isp < nsp; ++isp) neosp_i[isp+1] = 0;

    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        ++neosp_i[spost[ist][i]+1];
      }
    }

    neosp_i[0] = 0;
    for (int isp = 0; isp < nsp; ++isp) neosp_i[isp+1] += neosp_i[isp];
    assert(neosp_i[nsp] == nst*3);

    uint8 (*neosp_v)[2] = new uint8[neosp_i[nsp]][2];
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int nos = neosp_i[spost[ist][i]]++;
        neosp_v[nos][0] = spost[ist][(i+1)%3];
        neosp_v[nos][1] = (uint8(i)<<EDGE_OFFSET) | ist;
      }
    }

    // rewind...
    for (int isp = nsp-1; isp > 0; --isp) neosp_i[isp] = neosp_i[isp-1];
    neosp_i[0] = 0;

    for (int ist = 0; ist < nst; ++ist) FOR_I3 teost_data[ist][i] = 0; // no bits indicates an open edge

    // used inside the loop...
    vector<uint8> te_10_vec; // this is the vector of tri/edge nbrs where the matching edge is flipped
    // (i.e. 10) relative to our 01 edge -- the usual case
    vector<uint8> te_01_vec; // this is the vector of tri/edge nbrs where the matching edge is the same
    // (i.e. 01) as our 01 edge -- the tri nbr is probably misaligned in terms of its normal

    for (int isp0 = 0; isp0 < nsp; ++isp0) {
      for (int nos0 = neosp_i[isp0]; nos0 != neosp_i[isp0+1]; ++nos0) {
        // recall the neosp_v[nos0][0] is isp1 for this edge, but we don't need it yet.
        // the teost_data[ist0][i0] is used to check if we have already visited this teost...
        const int ist0 = neosp_v[nos0][1] & MASK_VALUE; const int i0 = (neosp_v[nos0][1]>>EDGE_OFFSET) & MASK_EDGE;
        if (teost_data[ist0][i0] == 0) {
          // we have never visited this teost before, so build vectors of ALL nbrs...
          assert(te_10_vec.empty());
          assert(te_01_vec.empty());
          // we have found a teost that has not been visited...
          // start by continuing through the loop looking for any misaligned matches...
          const int isp1 = neosp_v[nos0][0];
          for (int nos0b = nos0+1; nos0b != neosp_i[isp0+1]; ++nos0b) {
            if (neosp_v[nos0b][0] == uint8(isp1)) {
              // this is an edge of another tri that is aligned the same as our edge...
              te_01_vec.push_back(neosp_v[nos0b][1]);
            }
          }
          // and the properly opposite-aligned ones...
          for (int nos1 = neosp_i[isp1]; nos1 != neosp_i[isp1+1]; ++nos1) {
            if (neosp_v[nos1][0] == uint8(isp0)) {
              // we got a match with proper orientation....
              te_10_vec.push_back(neosp_v[nos1][1]);
            }
          }
          // we have all nbrs, so now process the vecs...
          if ((te_10_vec.size() == 1) && te_01_vec.empty()) {
            // this is the usual case...
            teost_data[ist0][i0] = (uint8(1)<<BITS_OFFSET) | te_10_vec[0];
            // and set the edge in te_10_vec[0] to point back to us...
            const int ist1 = te_10_vec[0] & MASK_VALUE; const int i1 = (te_10_vec[0]>>EDGE_OFFSET) & MASK_EDGE;
            assert(teost_data[ist1][i1] == 0);
            teost_data[ist1][i1] = (uint8(1)<<BITS_OFFSET) | neosp_v[nos0][1];
            te_10_vec.clear();
          }
          else if (te_10_vec.empty() && (te_01_vec.size() == 1)) {
            // this is the case of a flipped match...
            teost_data[ist0][i0] = (uint8(2)<<BITS_OFFSET) | te_01_vec[0];
            const int ist1 = te_01_vec[0] & MASK_VALUE; const int i1 = (te_01_vec[0]>>EDGE_OFFSET) & MASK_EDGE;
            assert(teost_data[ist1][i1] == 0);
            teost_data[ist1][i1] = (uint8(2)<<BITS_OFFSET) | neosp_v[nos0][1];
            te_01_vec.clear();
          }
          else if (!te_10_vec.empty() || !te_01_vec.empty()) {
            // here we have more than 2 tris meeting at the same edge -- a multiedge!
            teost_data[ist0][i0] = (uint8(3)<<BITS_OFFSET) | noome_vec.size(); // store the index of the multiedge, with 00 in 60,61 bits to indicate alignment
            teome_v_vec.push_back( neosp_v[nos0][1] );  // orientation in 62/63 is zero
            for (int ii = 0,ii_end=te_01_vec.size(); ii < ii_end; ++ii) {
              // each 01 gets the same conection...
              const int ist1 = te_01_vec[ii] & MASK_VALUE;
              const int i1 = (te_01_vec[ii]>>EDGE_OFFSET) & MASK_EDGE;
              assert(teost_data[ist1][i1] == 0);
              teost_data[ist1][i1] = (uint8(3)<<BITS_OFFSET) | noome_vec.size();
              teome_v_vec.push_back( te_01_vec[ii] ); // orientation in 62/63 is zero
            }
            for (int ii = 0,ii_end=te_10_vec.size(); ii < ii_end; ++ii) {
              const int ist1 = te_10_vec[ii] & MASK_VALUE;
              const int i1 = (te_10_vec[ii]>>EDGE_OFFSET) & MASK_EDGE;
              assert(teost_data[ist1][i1] == 0);
              teost_data[ist1][i1] = (uint8(3)<<BITS_OFFSET) | (uint8(1)<<EDGE_OFFSET) | noome_vec.size(); // orientation flipped wrt multiedge
              teome_v_vec.push_back( (uint8(1)<<BITS_OFFSET) | te_10_vec[ii] ); // orientation in 62/63 is 1 (i.e. fliped wrt multiedge)
            }
            noome_vec.push_back(pair<int,int>(isp0,isp1));
            teome_i_vec.push_back(teome_v_vec.size());
            // and clear...
            te_10_vec.clear();
            te_01_vec.clear();
          }
        }
      }
    }

    DELETE(neosp_v);
    DELETE(neosp_i);

    // check multiedges...

    for (int ime = 0,nme=noome_vec.size(); ime < nme; ++ime) {
      //cout << "ime: " << ime << " nodes: " << noome_vec[ime].first << " " << noome_vec[ime].second << " number of tris: " << teome_i_vec[ime+1]-teome_i_vec[ime] << endl;
      // check teost's are set properly...
      for (int tom = teome_i_vec[ime]; tom != teome_i_vec[ime+1]; ++tom) {
        // teome_v_vec[tom] can be unpacked as follows...
        const int orient = (teome_v_vec[tom]>>BITS_OFFSET); assert((orient == 0)||(orient == 1));
        const int i = (teome_v_vec[tom]>>EDGE_OFFSET) & MASK_EDGE; assert((i >= 0)&&(i < 3));
        const int ist = teome_v_vec[tom] & MASK_VALUE;
        assert(teost_data[ist][i] == ( (uint8(3)<<BITS_OFFSET) | (uint8(orient)<<EDGE_OFFSET) | ime ));
        if (orient == 0) {
          assert(spost[ist][i] == noome_vec[ime].first);
          assert(spost[ist][(i+1)%3] == noome_vec[ime].second);
        }
        else {
          assert(orient == 1);
          assert(spost[ist][i] == noome_vec[ime].second);
          assert(spost[ist][(i+1)%3] == noome_vec[ime].first);
        }
      }
    }

    // add imprinted edges
    // COUT2("build nooie.size: " << nooie_set.size());
    // if (nooie_set.size()) {
    //   int8 count = 0;
    //   set<edgeNodes>::iterator it;
    //   for (int ist = 0; ist < nst; ++ist) {
    //     FOR_I3 {
    //       it = nooie_set.find(edgeNodes(xsp[spost[ist][i]],xsp[spost[ist][(i+1)%3]]));
    //       if (it != nooie_set.end()) {
    //         teost_data[ist][i] |= (uint8(1)<<IMPRINT_OFFSET);
    //         ++count;
    //       }
    //     }
    //   }
    //   COUT2(" > imprinted half-edges: " << count);
    // }

    // report teost stats

    int8 counts[4] = { 0, 0, 0, 0 };
    for (int ist = 0; ist < nst; ++ist) {
      FOR_I3 {
        const int bits = (teost_data[ist][i]>>BITS_OFFSET); assert((bits >= 0)&&(bits < 4));
        ++counts[bits];
      }
    }

    COUT1(" > half-edge summary: open: " << counts[0] << ", connected: " << counts[1] << ", connected but flipped: " << counts[2] << ", multi: " << counts[3]);

    b_non_manifold[0] = (counts[0]) ? true:false;
    b_non_manifold[1] = (counts[2]) ? true:false;
    b_non_manifold[2] = (counts[3]) ? true:false;

    return true;
  }

  // Accessors to edge information

  int isManifold() const {
    if (teost_data == NULL) return -1;
    else {
      return (b_non_manifold[0] || b_non_manifold[1] || b_non_manifold[2]) ? 0:1;
    }
  }

  bool hasOpen() const {
    return b_non_manifold[0];
  }

  bool hasMisaligned() const {
    return b_non_manifold[1];
  }

  bool hasMulti() const {
    return b_non_manifold[2];
  }

  inline int getBits(const int ist,const int i) const {
    return int(teost_data[ist][i]>>BITS_OFFSET);
  }

  inline void checkValidInput(const int ist,const int i) const {
    assert(teost_data);
    assert((ist >= 0)&&(ist < nst));
    assert((i >= 0)&&(i < 3));
  }

  bool getTriNbrData(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);

    if ((bits == 1)||(bits == 2)) {
      ist_nbr = teost_data[ist][i] & MASK_VALUE; assert((ist_nbr >= 0)&&(ist_nbr < nst));
      i_nbr   = (teost_data[ist][i]>>EDGE_OFFSET) & MASK_EDGE; assert((i_nbr >= 0)&&(i_nbr < 3));
      orient_nbr = bits-1;
      return true;
    }
    return false;
  }

  int getTriNbrDataFull(int& ist_nbr,int& i_nbr,int& orient_nbr,const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);

    if ((bits == 1)||(bits == 2)) {
      ist_nbr = teost_data[ist][i] & MASK_VALUE; assert((ist_nbr >= 0)&&(ist_nbr < nst));
      i_nbr   = (teost_data[ist][i]>>EDGE_OFFSET) & MASK_EDGE; assert((i_nbr >= 0)&&(i_nbr < 3));
      orient_nbr = bits-1;
    }
    else if (bits == 3) {
      ist_nbr = teost_data[ist][i] & MASK_VALUE; assert((ist_nbr >= 0)&&(ist_nbr < int(noome_vec.size())));
      orient_nbr = (teost_data[ist][i]>>EDGE_OFFSET) & MASK_ME_ORIENT; assert((orient_nbr >= 0)&&(orient_nbr < 2));
    }
    assert((bits>=0) && (bits <= 3));
    return bits;
  }

  void setTriNbrData(const int ist_nbr,const int i_nbr,const int orient_nbr,const int ist,const int i) {
    checkValidInput(ist,i);

    assert((ist_nbr >= 0)&&(ist_nbr < nst));
    assert((i_nbr >= 0)&&(i_nbr < 3));
    assert((orient_nbr == 0)||(orient_nbr == 1));
    teost_data[ist][i] = (uint8(orient_nbr+1)<<BITS_OFFSET) | (uint8(i_nbr)<<EDGE_OFFSET) | ist_nbr;
  }

  void setTriNbrOpen(const int ist,const int i) {
    checkValidInput(ist,i);
    teost_data[ist][i] = 0;  // unindexed open edge
  }

  void setTriNbrMulti(const int ime,const int orient_ime,const int ist,const int i) {
    checkValidInput(ist,i);
    assert((orient_ime >= 0)&&(orient_ime < 2));
    assert((ime >= 0)&&(ime < int(noome_vec.size())));

    teost_data[ist][i] = (uint8(3)<<BITS_OFFSET) | (uint8(orient_ime)<<EDGE_OFFSET) | ime;
  }

  void setTriNbrToSameAs(const int ist_ref,const int i_ref,const int ist,const int i,const bool b_flip_orient) {
    const int flip_orient = (b_flip_orient) ? 1:0;

    checkValidInput(ist,i);
    checkValidInput(ist_ref,i_ref);
    const int bits = getBits(ist_ref,i_ref);
    if (bits == 0) {
      teost_data[ist][i] = teost_data[ist_ref][i_ref];  // open edges don't care about orientation...
    }
    else if ((bits == 1)||(bits == 2)) {
      const int ist_nbr = teost_data[ist_ref][i_ref] & MASK_VALUE; assert((ist_nbr >= 0)&&(ist_nbr < nst));
      const int i_nbr   = (teost_data[ist_ref][i_ref]>>EDGE_OFFSET) & MASK_EDGE; assert((i_nbr >= 0)&&(i_nbr < 3));
      const int orient_nbr = (bits-1+flip_orient)%2;
      teost_data[ist][i] = (uint8(orient_nbr+1)<<BITS_OFFSET) | (uint8(i_nbr)<<EDGE_OFFSET) | ist_nbr;
    }
    else {
      assert(bits == 3);
      const int ime = teost_data[ist_ref][i_ref] & MASK_VALUE; assert((ime >= 0)&&(ime < int(noome_vec.size())));
      int orient_ime = (teost_data[ist_ref][i_ref]>>EDGE_OFFSET) & MASK_ME_ORIENT; assert((orient_ime >= 0)&&(orient_ime < 2));
      orient_ime = (orient_ime+flip_orient)%2;
      teost_data[ist][i] = (uint8(bits)<<BITS_OFFSET) | (uint8(orient_ime)<<EDGE_OFFSET) | ime;
    }
  }

  void getTriNbrs(int (&ist_nbr)[3],const int ist) const {
    checkValidInput(ist,0);
    FOR_I3 {
      const int bits = getBits(ist,i);

      if ((bits == 1)||(bits == 2)) {
        ist_nbr[i] = teost_data[ist][i] & MASK_VALUE; assert((ist_nbr[i] >= 0)&&(ist_nbr[i] < nst));
      }
      else ist_nbr[i] = -1;  // for now treat me as an open edge...
    }
  }

  bool getAlignedTriNbr(int& ist_nbr,const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);
    if (bits == 1) {
      ist_nbr = teost_data[ist][i] & MASK_VALUE; assert((ist_nbr >= 0)&&(ist_nbr < nst));
      return true;
    }
    return false;
  }

  bool isNbrAligned(const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);
    if (bits == 1) return true;
    else return false;
  }

  bool isEdgeOpen(const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);

    return bits == 0;
  }

  bool isEdgeMulti(const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);

    return (bits == 3) ? true:false;
  }

  bool isEdgeMulti(int& ime,int& orient_ime,const int ist,const int i) const {
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);

    if (bits == 3) {
      ime = teost_data[ist][i] & MASK_VALUE; assert((ime >= 0)&&(ime < int(noome_vec.size())));
      orient_ime = (teost_data[ist][i]>>EDGE_OFFSET) & MASK_2BITS; assert((orient_ime >= 0)&&(orient_ime < 2));
      return true;
    }
    else return false;
  }

  void getMultiNbrs(vector<pair<int,int> >& ist_nbrs, const int ist,const int i) const {
    const int bits = getBits(ist,i);

    if (bits == 3) {
      const int ime = teost_data[ist][i] & MASK_VALUE; assert((ime >= 0)&&(ime < int(noome_vec.size())));
      for (int tom = teome_i_vec[ime]; tom != teome_i_vec[ime+1]; ++tom) {
        const int ist_nbr = teome_v_vec[tom] & MASK_VALUE;
        const int orient = (teome_v_vec[tom]>>BITS_OFFSET); assert((orient == 0)||(orient == 1));
        ist_nbrs.push_back(pair<int,int> (ist_nbr,orient));
      }
    }
  }

  void getMultiNbrsEdgeIndices(vector<pair<int,int> >& ist_nbrs, const int ist,const int i) const {
    const int bits = getBits(ist,i);
    if (bits == 3) {
      const int ime = teost_data[ist][i] & MASK_VALUE; assert((ime >= 0)&&(ime < int(noome_vec.size())));
      for (int tom = teome_i_vec[ime]; tom != teome_i_vec[ime+1]; ++tom) {
        const int ist_nbr = teome_v_vec[tom] & MASK_VALUE;
        const int i_nbr =(teome_v_vec[tom]>>EDGE_OFFSET) & MASK_EDGE; assert((i_nbr >= 0)&&(i_nbr < 3));
        ist_nbrs.push_back(pair<int,int> (ist_nbr,i_nbr));
      }
    }
  }

  void getMultiFaces(vector<int>& faces, int ime) const {
    for (int tom = teome_i_vec[ime]; tom != teome_i_vec[ime+1]; ++tom) {
      const int ist_nbr = teome_v_vec[tom] & MASK_VALUE;
      faces.push_back(ist_nbr);
    }
  }

  void setOpenEdgeGroup(const int igr,const int ist,const int i) {
    assert(igr >= 0);
    checkValidInput(ist,i);
    const int bits = getBits(ist,i);
    assert(bits == 0);
    teost_data[ist][i] = igr;
  }

  bool getOpenEdgeGroup(int& igr,const int ist,const int i) const {

    checkValidInput(ist,i);
    const int bits = getBits(ist,i);
    if (bits == 0) {
      igr = int(teost_data[ist][i]);
      return true;
    }
    return false;
  }

  bool isEdgeSelectedBoundary(const int ist, const int i,const IntFlag& st_flag) const {
    // any edges who don't have neighbor with st_flag == 1 are boundaries
    if (st_flag[ist] == 0) return false;

    int ist_nbr,i_nbr,orient_nbr;
    const bool valid_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
    if (!valid_nbr) {
      int ime,orient_ime;
      if(isEdgeMulti(ime,orient_ime,ist,i)) {
        vector<int> ist_nbrs;
        getMultiFaces(ist_nbrs,ime);
        for (vector<int>::iterator ist_nbr=ist_nbrs.begin(); ist_nbr!=ist_nbrs.end(); ++ist_nbr) {
          if (st_flag[*ist_nbr] == 0) return true;  // if any multi-nbr is not flagged, consider it a boundary...?
        }
        return false;  // multi with all flagged
      }
      return true;  // open should be only choice
    }
    else {
      if (st_flag[ist_nbr] == 0) return true;
      else return false;
    }
  }

  bool isEdgeZoneBoundary(const int ist, const int i,const int *znost) const {
    // open edges of zone currently treated differently from zone boundary. These are edges between two zones.
    int ist_nbr,i_nbr,orient_nbr;
    const bool valid_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
    if (!valid_nbr) return false;
    else {
      if (znost[ist] == znost[ist_nbr]) return false;
      else return true;
    }
  }

  bool isEdgeZoneBoundary(const int ist, const int i,const IntFlag& znost) const {
  // open edges of zone currently treated differently from zone boundary. These are edges between two zones.
  int ist_nbr,i_nbr,orient_nbr;
  const bool valid_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
  if (!valid_nbr) {
    int ime,orient_ime;
    if(isEdgeMulti(ime,orient_ime,ist,i)) {
      vector<int> ist_nbrs;
      getMultiFaces(ist_nbrs,ime);
      for (vector<int>::iterator ist_nbr=ist_nbrs.begin(); ist_nbr!=ist_nbrs.end(); ++ist_nbr) {
        if (znost[ist] != znost[*ist_nbr]) return true;  // if any multi-nbr is not my zone, consider it a boundary...
      }
    }

    return false;  // either an multi with all same zone OR an open edge...is this considered the zone boundary too?
    }
    else {
      if (znost[ist] == znost[ist_nbr]) return false;
      else return true;
    }
  }

  bool getZoneBoundaryPair(pair<int,int>& zonePair,const int ist, const int i,const int *znost) const {
    // open edges of zone currently treated differently from zone boundary. These are edges between two zones.
    int ist_nbr,i_nbr,orient_nbr;
    const bool valid_nbr = getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i);
    if (!valid_nbr) return false;
    else {
      if (znost[ist] == znost[ist_nbr]) return false;

      zonePair.first = min(znost[ist],znost[ist_nbr]);
      zonePair.second = max(znost[ist],znost[ist_nbr]);
      return true;
    }
  }

  // void setImprintedHalfEdge(const double isp0[3], const double isp1[3]) {
  //   nooie_set.insert(edgeNodes (DoubleVertex(isp0),DoubleVertex(isp1)));
  // }
  //
  // bool isEdgeImprinted(const int ist,const int i) const {
  //   checkValidInput(ist,i);
  //   return (int( (teost_data[ist][i]>>IMPRINT_OFFSET) & MASK_IMPRINT ) == 1) ? true:false;
  // }
  //
  // void clearImprintedHalfEdges() {
  //   nooie_set.clear();
  // }

};

#undef BITS_OFFSET
#undef EDGE_OFFSET
#undef MASK_VALUE
#undef MASK_BITS

#undef MASK_EDGE
#undef MASK_ME_ORIENT

#endif
