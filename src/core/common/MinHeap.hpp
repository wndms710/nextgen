#ifndef _MIN_HEAP_HPP_
#define _MIN_HEAP_HPP_

//#include "CTI.hpp"
//using namespace CTI;

// maybe template this in the future
class MinHeap {

private:

  vector<pair<double,int> > heap; // vector data
  int nsp;
  int * index_sp; // back pointers to grid

  // return left child
  int leftChild(int parent) {
    int l = 2 * parent + 1;
    if (l < int(heap.size())) return l;
    else return -1;
  }

  // return right child
  int rightChild(int parent) {
    int r = 2 * parent + 2;
    if (r < int(heap.size())) return r;
    else return -1;
  }

  // return parent
  int parent(int child) {
    int p = (child - 1)/2;
    if (child == 0) return -1;
    else return p;
  }

  // maintain heap structure from the bottom up
  void upHeap(int in) {
    if (in >= 0 && parent(in) >= 0 && heap[parent(in)] > heap[in]) {
      pair<double,int> temp = heap[in];
      heap[in] = heap[parent(in)];
      heap[parent(in)] = temp;
      // update index_sp
      if (index_sp) {
        index_sp[heap[in].second] = in;
        index_sp[heap[parent(in)].second] = parent(in);
      }
      // next level
      upHeap(parent(in));
    }
  }

  // maintain heap structure from top down
  void downHeap(int in) {
    int child = leftChild(in);
    int child1 = rightChild(in);
    if (child >= 0 && child1 >= 0 && heap[child] > heap[child1]) {
      child = child1;
    }
    if (child > 0 && heap[in] > heap[child]) {
      pair<double,int> temp = heap[in];
      heap[in] = heap[child];
      heap[child] = temp;
      // update index_sp
      if (index_sp) {
        index_sp[heap[in].second] = in;
        index_sp[heap[child].second] = child;
      }
      // next level
      downHeap(child);
    }
  }

public:

  // constructor
  MinHeap() {
    index_sp = NULL;
  }

  // constructor
  MinHeap(int nsp) {
    index_sp = new int[nsp];
    for (int isp = 0; isp < nsp; ++isp) index_sp[isp] = -1;
    this->nsp = nsp;
  }

  // destructor
  ~MinHeap() {
    heap.clear();
    DELETE(index_sp);
  }

  // heap size
  int size() {
    return heap.size();
  }

  // insert element into heap
  void insert(pair<double,int> element)
  {
    if (index_sp) index_sp[element.second] = heap.size();
    heap.push_back(element);
    upHeap(heap.size()-1);
  }

  // delete minimum element from heap
  void deleteMin()
  {
    if (heap.size() == 0) return;
    if (index_sp) index_sp[heap[heap.size()-1].second] = -1;
    heap[0] = heap[heap.size()-1];
    heap.pop_back();
    downHeap(0);
  }

  // extract minimum element from heap
  pair<double,int> extractMin() {
    assert(heap.size() > 0);
    return heap.front();
  }

  // display heap for debugging
  void displayHeap() {
    vector<pair<double,int> >::iterator iter = heap.begin();
    while (iter != heap.end()) {
      cout << "(" << iter->first << "," << iter->second << ") ";
      ++iter;
    }
    cout << endl;
  }

  // display indices for debugging
  void displayIndices() {
    if (index_sp) {
      for (int isp = 0; isp < nsp; ++isp) {
        if (index_sp[isp] >= 0) {
           cout <<"(" << index_sp[isp] << "," << isp << ") ";
        }
      }
      cout << endl;
    }
  }

  // takes min of element in double and input double and then updates heap
  void updateDouble(int in, double dist) {
    heap[in].first = min(heap[in].first,dist);
    upHeap(in);
  }

  // returns index to element whose pair.second = isp
  int index(int isp) {
    return index_sp[isp];
  }
};

#endif
