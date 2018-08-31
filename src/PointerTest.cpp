#include<stdio.h>
#include<ostream>
using namespace std;
void changePointer(int* a);
int main(){
  int* a = new int[2];
  a[0] = 1;
  a[1] = 2;
  changePointer(a);
  //cout << a[0] << "," << a[1] << endl;
  printf("%d,%d\n", a[0], a[1]);
  return 0;
}
void changePointer(int* a){
  //a[0] ++;
  //a[1] ++;
  //
  *a += 1;
  *(a+1) += 1;
}
