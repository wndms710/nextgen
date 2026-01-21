#ifndef _NONMANIFOLD_TRI_FLAG_HPP_
#define _NONMANIFOLD_TRI_FLAG_HPP_

#include "CTI.hpp"
#include "IntFlag.hpp"
using namespace CTI;

enum NM_FLAG {
  NM_MULTI_EDGE_ADJACENT = (1<<0),
  NM_COLLAPSED = (1<<1),
  NM_LINEAR = (1<<2),
  NM_SHARED_EDGES = (1<<3),
  NM_INTERSECT = (1<<4),
  NM_OVERLAP = (1<<5),
  NM_IMPRINT = (1<<6),
  NM_VISITED = (1<<7)
};

class NonManifoldTriFlag {
// bit-setting based wrapper around an IntFlag; bits indicate different states
private:
  IntFlag my_st_flag;
  IntFlag n_counts;
  double dist_tol;
  double angle_tol_degrees;
  bool check_self_intersections;
  string output_format;
  int n_seed_samples;
  int nbr_layers;

public:
  NonManifoldTriFlag() {
    n_counts.setLength(7);
    resetCounts();
    dist_tol = 1.0e-10;
    angle_tol_degrees = 0.5;
    check_self_intersections = false;
    output_format = "none";
    n_seed_samples = nbr_layers = 0;
  };

  NonManifoldTriFlag(const int size) {
    my_st_flag.resize(size);
    n_counts.setLength(7);
    resetCounts();
    dist_tol = 1.0e-10;
    angle_tol_degrees = 0.5;
    check_self_intersections = false;
    output_format = "none";
    n_seed_samples = nbr_layers = 0;
  };

  ~NonManifoldTriFlag() {
    clear();
  }
  void clear() {
    my_st_flag.clear();
    resetCounts();
  }
  void resetCounts() {
    n_counts.setAll(0);
  }
  bool isNull() const {
    return my_st_flag.isNull();
  }
  void setDistTol(const double _dist_tol) {
    dist_tol = _dist_tol;
  }
  double distTol() const {
    return dist_tol;
  }
  void setAngleTol(const double _angle_tol) {
    angle_tol_degrees = _angle_tol;
  }
  double angleTol() const {
    return angle_tol_degrees;
  }
  void setCheckSelfIntersections(const bool _check) {
    check_self_intersections = _check;
  }
  bool checkSelfIntersections() const {
    return check_self_intersections;
  }
  void setFormat(string _format) {
    transform(_format.begin(), _format.end(), _format.begin(), ::tolower);
    if (_format == "tecplot" || _format == "sbin" || _format == "none") output_format = _format;
    else CWARN("unrecognized output format \"" << _format << "\"; producing no output");
  }
  string getFormat() {
    return output_format;
  }
  void setNSeedSamples(const int n) {
    n_seed_samples = n;
  }
  int nSeedSamples() const {
    return n_seed_samples;
  }
  void setNbrLayers(const int n) {
    nbr_layers = n;
  }
  int nNbrLayers() const {
    return nbr_layers;
  }

  void dump() {
    if (n_counts[0]) COUT1(" > non-manifold edge adjacent tris: " << n_counts[0]);
    if (n_counts[1]) COUT1(" > collapsed tris: " << n_counts[1]);
    if (n_counts[2]) COUT1(" > linear tris: " << n_counts[2]);
    if (n_counts[3]) COUT1(" > multi-edge neighbor tris: " << n_counts[3]);
    if (n_counts[4]) COUT1(" > intersecting tris: " << n_counts[4]);
    if (n_counts[5]) COUT1(" > overlapping tris: " << n_counts[5]);
    if (n_counts[6]) COUT1(" > imprinting tris: " << n_counts[6]);

    if (n_counts.count() == 0) COUT1(" > no issues currently registered; (re)run RUN_DIAGNOSTICS to check for potential issues.");
  }
  void getCounts() {
    // allows intialization of counts from an arbitrary st_flag
    resetCounts();
    for (int ist=0; ist < my_st_flag.size(); ++ist) {
      if (my_st_flag[ist]&NM_MULTI_EDGE_ADJACENT) ++n_counts[0];
      if (my_st_flag[ist]&NM_COLLAPSED) ++n_counts[1];
      if (my_st_flag[ist]&NM_LINEAR) ++n_counts[2];
      if (my_st_flag[ist]&NM_SHARED_EDGES) ++n_counts[3];
      if (my_st_flag[ist]&NM_INTERSECT) ++n_counts[4];
      if (my_st_flag[ist]&NM_OVERLAP) ++n_counts[5];
      if (my_st_flag[ist]&NM_IMPRINT) ++n_counts[6];
    }
  }

  // accessors and setters
  inline int sumCounts() const {
    return n_counts.count();
  }
  inline void setAll(const int val) {
    my_st_flag.setAll(val);
    getCounts();
  }
  inline void resize(const int n) {
    my_st_flag.resize(n);
  }
  inline bool isClean (const int ist) const {
    return ((my_st_flag[ist]&(~NM_VISITED)) == 0) ? true:false;  // uses any
  }

  inline void setClean(const int ist) {
    my_st_flag[ist] = 0;
  }

  inline bool hasBits(const int ist, const int nme_bits) {
    return (my_st_flag[ist] & nme_bits);
  }

  inline void setNmeAdjacent(const int ist) {
    my_st_flag[ist] |= NM_MULTI_EDGE_ADJACENT;
  }
  inline void unSetNmeAdjacent(const int ist) {
    my_st_flag[ist] &= ~NM_MULTI_EDGE_ADJACENT;
  }
  inline bool isNmeAdjacent(const int ist) const {
    return (my_st_flag[ist] & NM_MULTI_EDGE_ADJACENT);
  }
  inline int countNmeAdjacent() {
    return n_counts[0];
  }
  inline void setCollapsed(const int ist) {
    my_st_flag[ist] |= NM_COLLAPSED;
  }
  inline void unsetCollapsed(const int ist) {
    my_st_flag[ist] &= ~NM_COLLAPSED;
  }
  inline bool isCollapsed(const int ist) const {
    return (my_st_flag[ist] & NM_COLLAPSED);
  }
  inline int countCollapsed() {
    return n_counts[1];
  }
  inline void setLinear(const int ist) {
    my_st_flag[ist] |= NM_LINEAR;
  }
  inline void unsetLinear(const int ist) {
    my_st_flag[ist] &= ~NM_LINEAR;
  }
  inline bool isLinear(const int ist) const {
    return (my_st_flag[ist] & NM_LINEAR);
  }
  inline int countLinear() {
    return n_counts[2];
  }
  inline void setMultiNeighbor(const int ist) {
    my_st_flag[ist] |= NM_SHARED_EDGES;
  }
  inline void unsetMultiNeighbor(const int ist) {
    my_st_flag[ist] &= ~NM_SHARED_EDGES;
  }
  inline bool isMultiNeighbor(const int ist) const {
    return (my_st_flag[ist] & NM_SHARED_EDGES);
  }
  inline int countMultiNeighbor() {
    return n_counts[3];
  }
  inline void setIntersecting(const int ist) {
    my_st_flag[ist] |= NM_INTERSECT;
  }
  inline void unsetIntersecting(const int ist) {
    my_st_flag[ist] &= ~NM_INTERSECT;
  }
  inline bool isIntersecting(const int ist) const {
    return (my_st_flag[ist] & NM_INTERSECT);
  }
  inline int countIntersecting() {
    return n_counts[4];
  }
  inline void setOverlapping(const int ist) {
    my_st_flag[ist] |= NM_OVERLAP;
  }
  inline void unsetOverlapping(const int ist) {
    my_st_flag[ist] &= ~NM_OVERLAP;
  }
  inline bool isOverlapping(const int ist) const {
    return (my_st_flag[ist] & NM_OVERLAP);
  }
  inline int countOverlapping() {
    return n_counts[5];
  }
  inline void setImprinting(const int ist) {
    my_st_flag[ist] |= NM_IMPRINT;
  }
  inline void unsetImprinting(const int ist) {
    my_st_flag[ist] &= ~NM_IMPRINT;
  }
  inline bool isImprinting(const int ist) const {
    return (my_st_flag[ist] & NM_IMPRINT);
  }
  inline int countImprinting() {
    return n_counts[6];
  }
  inline void setVisited(const int ist) {
    my_st_flag[ist] |= NM_VISITED;
  }
  inline bool isVisited(const int ist) const {
    return (my_st_flag[ist] & NM_VISITED);
  }
  inline void unsetVisited(const int ist) {
    my_st_flag[ist] &= ~NM_VISITED;
  }
};

#endif
