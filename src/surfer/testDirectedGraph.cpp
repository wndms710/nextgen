#include "DirectedGraph.hpp"

// Driver program to test above functions
int main() {

  // Create a graph given in the above diagram
  DirectedGraph g(9);

  // made 1 indexed to test out having different internal and external edge and vertex indices
  g.addEdge(8,9,1);
  g.addEdge(9,8,2);

  g.addEdge(2,9,3);
  g.addEdge(1,8,4);
  g.addEdge(2,7,5);

  g.addEdge(1,2,6);
  g.addEdge(3,1,7);
  g.addEdge(3,2,8);
  g.addEdge(2,3,9);
  g.addEdge(3,4,10);
  g.addEdge(3,6,11);
  g.addEdge(6,4,12);
  g.addEdge(4,5,13);
  g.addEdge(5,2,14);

  g.addEdge(1,5,15);

  //g.runStronglyConnectedComponentFinder(); // called by runElementaryCircuitFinder
  //g.runElementaryCircuitFinder(); // called by runChordlessCircuitFinder
  g.runChordlessCircuitFinder();

  cout << "Coordless Circuit edges: " << endl;
  const int ncc = g.getNumberOfChordlessCircuits();
  for (int icc = 0; icc < ncc; ++icc) {
    int ned = g.getNumberOfChordlessCircuitEdges(icc);
    cout << icc <<": ";
    for (int ied = 0; ied < ned; ++ied) {
      cout << g.getChordlessCircuitEdge(icc,ied) << " ";
    }
    cout << endl;
  }
  return 0;
}

