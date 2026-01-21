#include "DirectedGraph.hpp"

// A recursive function to perform DFS starting from v, fills scc
void DirectedGraph::findSCCVertices(const int v) {

  // mark the current node as visited and index it
  visited[v] = true;
  scc[v] = nscc;

  // recur for all the vertices adjacent to this vertex
  for (int w = 0; w < nv; ++w) {
    if (adj[v][w] >= 0) {
      if (!visited[w]) {
        findSCCVertices(w);
      }
    }
  }
}

// fills v_stack with vertices (in increasing order of finishing
// times). The top element of stack has the maximum finishing time
void DirectedGraph::stackVertices(stack<int>& v_stack, const int v) {

  // mark the current node as visited
  visited[v] = true;

  // recur for all the vertices adjacent to this vertex
  for (int w = 0; w < nv; ++w) {
    if (adj[v][w] >= 0) {
      if (!visited[w]) stackVertices(v_stack,w);
    }
  }

  // all vertices reachable from v are processed by now, push v
  v_stack.push(v);
}

// the main function that finds and prints all strongly connected components
void DirectedGraph::runStronglyConnectedComponentFinder() {

  // make sure graph is fully built
  assert(nvl == nv);
  iv_map.clear();

  // allocate some arrays and init pod for scc
  nscc = 0;
  visited = new bool[nv];
  scc = new int[nv];
  stack<int> v_stack;

  // mark all the vertices as not visited (for first DFS)
  for (int v = 0; v < nv; ++v)
    visited[v] = false;

  // should be >= 0 at the end (error checking)
  for (int v = 0; v < nv; ++v)
    scc[v] = -1;

  // fill vertices in stack according to their finishing times
  for (int v = 0; v < nv; ++v)
    if (!visited[v]) stackVertices(v_stack,v);

  // transpose adjacency matrix
  transposeAdjacencyMatrix();

  // mark all the vertices as not visited (for second DFS)
  for (int v = 0; v < nv; ++v)
    visited[v] = false;

  // now process all vertices in order defined by v_stack
  while (!v_stack.empty()) {

    // pop a vertex from stack
    int v = v_stack.top();
    v_stack.pop();

    // build strongly connected component associated to popped vertex
    if (!visited[v]) {
      findSCCVertices(v);
      ++nscc;
    }
  }

  // transpose adjacency matrix
  transposeAdjacencyMatrix();

  // print index and scc index
  cout << "The following are strongly connected component indices graph (v scc)" << endl;
  for (int v = 0; v < nv; ++v) cout << v_index[v] << " " << scc[v] << endl;

}

void DirectedGraph::initSCCSubGraph(const int iscc) {
  int vl = 0;
  for (int v = 0; v < nv; ++v) {
    if (scc[v] == iscc) {
      int wl = 0;
      for (int w = 0; w < nv; ++w) {
        if (scc[w] == iscc) {
          adj_scc[vl][wl] = (adj[v][w] >= 0); // only need bool here
          ++wl;
        }
      }
      iv_scc[vl] = v;
      blocked_scc[vl] = false;
      b_scc[vl].clear();
      ++vl;
    }
  }
  nv_scc = vl;
}

void DirectedGraph::unblock(int u) {
  blocked_scc[u] = false;
  while (!b_scc[u].empty()) {
    int w = b_scc[u].front();
    b_scc[u].pop_front();
    if (blocked_scc[w])
      unblock(w);
  }
}

bool DirectedGraph::circuit(vector<int>& v_vec_scc,int v,const int s) {

  bool f = false;
  v_vec_scc.push_back(v);
  blocked_scc[v] = true;

  for (int w = 0; w < nv_scc; ++w) {
    if (adj_scc[v][w]) {
      if (w == s) {
        addElementaryCircuit(v_vec_scc);
        f = true;
      }
      else if (w > s && !blocked_scc[w]) {
        f = circuit(v_vec_scc,w,s);
      }
    }
  }

  if (f) {
    unblock(v);
  }
  else {
    for (int w = 0; w < nv_scc; ++w) {
      if (adj_scc[v][w]) {
        list<int>::iterator it = find(b_scc[w].begin(), b_scc[w].end(), v);
        if (it == b_scc[w].end()) {
          b_scc[w].push_back(v);
        }
      }
    }
  }

  v_vec_scc.pop_back();
  return f;
}

void DirectedGraph::addElementaryCircuit(const vector<int>& v_vec_scc) {
  stringstream ss;
  ec_vec_v_vec.push_back(vector<int>(v_vec_scc.size()));
  for (int v = 0,v_max=v_vec_scc.size(); v < v_max; ++v) {
    ec_vec_v_vec[nec][v] = iv_scc[v_vec_scc[v]]; // local so we know how to index w/o map
    ss << v_index[iv_scc[v_vec_scc[v]]] << " "; // put actual indices for output purposed
  }
  ec_names.push_back(string());
  ec_names[nec] = ss.str() + ss.str(); // doubled to easily account for periodicity
  ++nec;
}

void DirectedGraph::runElementaryCircuitFinder() {

  // make sure we have sccs
  if (!visited) runStronglyConnectedComponentFinder();

  // scc structures use nv (can change to max number of nodes in scc's)
  blocked_scc = visited; // just renamed for clarity
  b_scc = new list<int>[nv];
  iv_scc = new int[nv];
  adj_scc = new bool*[nv];
  for (int i = 0; i < nv; ++i)
    adj_scc[i] = new bool[nv];
  nec = 0;
  vector<int> v_vec_scc;

  for (int iscc = 0; iscc < nscc; ++iscc) {
    initSCCSubGraph(iscc);

    v_vec_scc.clear();
    int s = 0;

    while (s < nv_scc-1) {
      for (int i = s; i < nv_scc; ++i) {
        blocked_scc[i] = false;
        b_scc[i].clear();
      }
      circuit(v_vec_scc,s,s);
      ++s;
    }
  }
  cout << "The following are elementary circuits in given strongly connected components" << endl;
  for (int iec = 0; iec < nec; ++iec)
    cout << ec_names[iec].substr(0,ec_names[iec].length()/2) << endl;
}

void DirectedGraph::runChordlessCircuitFinder() {

  // make sure we have elementary circuits
  if (!adj_scc) runElementaryCircuitFinder();

  ncc = 0;
  ec_chordless = new bool[nec];

  for (int iec = 0; iec < nec; ++iec) {
    ec_chordless[iec] = true;
    if (ec_vec_v_vec[iec].size() > 2) {
      for (int iec2 = 0; iec2 < nec; ++iec2) {
        if (iec2 != iec && ec_vec_v_vec[iec2].size() > 2) {
          // if I contain any body then i have a chord
          for (int i = 0,limit=ec_names[iec2].length()/2; i < limit; i+=2) {
            size_t found = ec_names[iec].find(ec_names[iec2].substr(i,ec_names[iec2].length()/2));
            if (found != string::npos) {
              ec_chordless[iec] = false;
              break;
            }
            if (!ec_chordless[iec]) break;
          }
        }
      }
    }
    if (ec_chordless[iec]) ++ncc;
  }

  ecocc = new int[ncc];
  int icc = 0;
  for (int iec = 0; iec < nec; ++iec) {
    if (ec_chordless[iec]) {
      ecocc[icc++] = iec;
    }
  }
  assert(icc == ncc);

  cout << "The following are chordless circuits among the elementary circuits" << endl;
  for (int iec = 0; iec < nec; ++iec)
    cout << ec_chordless[iec] << endl;
}

int DirectedGraph::getChordlessCircuitEdge(const int icc, const int ied) {
  // ied correspond to first v on edge
  int v = ec_vec_v_vec[ecocc[icc]][ied];
  int ied2 = ied+1;
  if (ied2 == int(ec_vec_v_vec[ecocc[icc]].size())) ied2 = 0;
  int w = ec_vec_v_vec[ecocc[icc]][ied2];
  return adj[v][w];
}
