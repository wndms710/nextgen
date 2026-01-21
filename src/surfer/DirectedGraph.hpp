#ifndef _DIRECTED_GRAPH_HPP_
#define _DIRECTED_GRAPH_HPP_

#include "CTI.hpp"
#include <stack>
#include <sstream>
using namespace CTI;

class DirectedGraph {
private:

  int nv; // number of vertices
  int **adj; // adjacency matrix (w/ edge (open edge group) label)
  int *v_index; // actual vertex index
  map<const int,int> iv_map; // actual vertex index to local index (can be cleared once graph is built)
  int nvl; //should = nv after adding all the edges

  // data structures for find strongly connected components (scc)
  bool *visited; // has this vertex been visited
  int nscc; // number of scc's
  int *scc; // scc index

  // data structures for finding elementary circuits within a scc
  int nv_scc; // number of vertices for current scc
  int *iv_scc; // mapping back to graph vertices
  bool **adj_scc; // adjacency matrix for current scc
  bool *blocked_scc; // has this vertex been blocked
  list<int> *b_scc; // helps avoid fruitless searches

  // store elementary circuits
  int nec;
  vector<vector<int> > ec_vec_v_vec; // circuits' local vertices

  // data structures for finding chordless circuits
  vector<string> ec_names;
  bool *ec_chordless;
  int ncc;
  int *ecocc;

  // transpose adjacency matrix
  void transposeAdjacencyMatrix() {
    for (int v = 0; v < nv; ++v) {
      for (int w = 0; w < v; ++w) {
        const int b = adj[w][v];
        adj[w][v] = adj[v][w];
        adj[v][w] = b;
      }
    }
  }

  // fills v_stack with vertices (in increasing order of finishing
  // times). the top element of stack has the maximum finishing time
  void stackVertices(stack<int>& v_stack,const int v);

  // recursive DFS starting from v, throwing unvisited vertices into scc
  void findSCCVertices(const int v);

  void initSCCSubGraph(const int iscc);
  void unblock(int u);
  bool circuit(vector<int>& v_vec_scc,int v,const int s);
  void addElementaryCircuit(const vector<int>& v_vec_scc);

public:

  DirectedGraph(int nv) {
    this->nv = nv;
    adj = new int*[nv];
    for (int i = 0; i < nv; ++i)
      adj[i] = new int[nv];
    for (int i = 0; i < nv; ++i)
      for (int j = 0; j < nv; ++j)
        adj[i][j] = -1;
    nvl = 0;
    v_index = new int[nv];
    for (int i = 0; i < nv; ++i)
      v_index[i] = -1; // for checks

    // these may or may not get allocated (nullify here)
    visited = NULL;
    scc = NULL;
    iv_scc = NULL;
    adj_scc = NULL;
    blocked_scc = NULL;
    b_scc = NULL;
    ec_chordless = NULL;
    ecocc = NULL;
  }

  ~DirectedGraph() {
    for (int i = 0; i < nv; ++i)
      delete[] adj[i];
    delete[] adj;
    delete[] v_index;

    // these may not have bee allocated
    if (visited) delete[] visited;
    if (scc) delete[] scc;
    if (adj_scc) {
      for (int i = 0; i < nv; ++i)
        delete[] adj_scc[i];
      delete[] adj_scc;
    }
    if (b_scc) delete[] b_scc;
    if (iv_scc) delete[] iv_scc;
    if (ec_chordless) delete[] ec_chordless;
    if (ecocc) delete[] ecocc;
  }

  // add v->w to adjacency matrix
  void addEdge(const int v, const int w, const int e) {
    if (nv == 0) return;
    map<const int,int>::iterator it = iv_map.find(v);
    int vl = -1;
    if (it == iv_map.end()) {
      vl = nvl++;
      iv_map[v] = vl;
      assert(v_index[vl] == -1);
      v_index[vl] = v;
    }
    else {
      vl = it->second;
    }
    it = iv_map.find(w);
    int wl = -1;
    if (it == iv_map.end()) {
      wl = nvl++;
      iv_map[w] = wl;
      assert(v_index[wl] == -1);
      v_index[wl] = w;
    }
    else {
      wl = it->second;
    }

    assert(vl >= 0 && wl >= 0);
    adj[vl][wl] = e;
  }

  // calc strongly connected components
  void runStronglyConnectedComponentFinder();

  // run circuit finder
  void runElementaryCircuitFinder();

  // get chordless circuits among circuits
  void runChordlessCircuitFinder();

  // accessors
  int getNumberOfChordlessCircuits() { return ncc; }
  int getNumberOfChordlessCircuitEdges(const int icc) { return ec_vec_v_vec[ecocc[icc]].size(); }
  int getChordlessCircuitEdge(const int icc, const int ied);
};

#endif
