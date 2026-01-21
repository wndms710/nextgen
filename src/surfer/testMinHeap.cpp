#include "MinHeap.hpp"

int main() {
  MinHeap h;
  while (1) {
    cout << "------------------" << endl;
    cout << "Operations on Heap" << endl;
    cout << "------------------" << endl;
    cout << "1. Insert Element" << endl;
    cout << "2. Delete Minimum Element" << endl;
    cout << "3. Extract Minimum Element" << endl;
    cout << "4. Edit Element" << endl;
    cout << "5. Print Heap" << endl;
    cout << "6. Exit" << endl;
    int choice, in_int;
    double in_double;
    pair<double,int> element;
    cout << "Enter your choice: ";
    cin >> choice;
    switch(choice) {
    case 1:
      cout << "Enter the double to be inserted: ";
      cin >> in_double;
      cout << "Enter the corresponding int to be inserted: ";
      cin >> in_int;
      h.insert(pair<double,int>(in_double,in_int));
      break;
    case 2:
      h.deleteMin();
      break;
    case 3:
      cout << "Minimum Element: ";
      if (h.size() == 0)
        cout << "Heap is Empty"<<endl;
      else
        element = h.extractMin();
        cout << "Minimum Element: (" << element.first << "," << element.second << ") " << endl;
      break;
    case 4:
      cout << "Enter the integer element index: ";
      cin >> in_int;
      cout << "Enter the corresponding double to be updated: ";
      cin >> in_double;
      h.updateDouble(in_int,in_double);
      break;
    case 5:
      cout << "Displaying elements of Heap: ";
      h.displayHeap();
      break;
    default:
      exit(1);
    }
  }
  return 0;
}
