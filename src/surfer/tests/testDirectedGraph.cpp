#include "../DirectedGraph.hpp"
#include "../../core/catch1/catch.hpp"

TEST_CASE( "directed graph","[directed_graph]" ) {
  INFO(" > generating dummy directed graph");
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

  const int ncc = g.getNumberOfChordlessCircuits();
  CHECK(ncc == 4);

  vector<vector<int> > results(ncc);

  for (int icc = 0; icc < ncc; ++icc) {
    int ned = g.getNumberOfChordlessCircuitEdges(icc);
    for (int ied = 0; ied < ned; ++ied) {
      results[icc].push_back(g.getChordlessCircuitEdge(icc,ied));
    }
  }

  // stdout for debugging of test
  // for (int icc=0; icc<ncc; ++icc) {
  //   cout << icc << ": ";
  //   for(vector<int>::iterator it=results[icc].begin(); it!=results[icc].end(); ++it) cout << *it << " ";
  // }
  // cout << endl;

  vector<int> answer;

  answer.clear();
  answer.push_back(9);
  answer.push_back(8);
  CHECK(results[0] == answer);

  answer.clear();
  answer.push_back(9);
  answer.push_back(7);
  answer.push_back(6);
  CHECK(results[1] == answer);

  answer.clear();
  answer.push_back(9);
  answer.push_back(10);
  answer.push_back(13);
  answer.push_back(14);
  CHECK(results[2] == answer);

  answer.clear();
  answer.push_back(1);
  answer.push_back(2);
  CHECK(results[3] == answer);
}

