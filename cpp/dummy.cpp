#include <iostream>
#include <vector>

using namespace std;

vector<int> dummy(int n) {
  vector<int> a(n, 1);
  return a;
}

void dummy2(int* p, int n) {
  for (int i=0; i < n; i++) {
    p[i] = 1;
  }
}

int main() {
  int* k = new int[4];
  k[0] = 1;
  k[1] = 2;
  k[2] = 3;
  k[3] = 4;

  int* s = new int;
  
  for (int i = 0; i<5; i++) {
    ++(*s);
  }
  vector<int> d = dummy(4);
  d = dummy(6);
  // for (int j=0; j < 4; j++){
  //   vector<int> h(k, k+4);
  //   for (int i=0; i < 4; i++) {
  //     cout << h[i] << endl;
  //   }
  // }
  int h[] = {3, 4};
  for (int i =0; i < 2; i++) {
    cout << h[i] << "a" << endl;
    int* da = new int[h[i]];
    dummy2(da, h[i]);
    for (int j=0; j < h[i]; j++) {
      cout << da[j] << endl;
    }
    cout << "-----\n";
  }


  return 0;
}