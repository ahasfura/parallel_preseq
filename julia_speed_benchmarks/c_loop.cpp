#include <iostream>
#include <cstdio>
#include <ctime>

int main(){
  clock_t start;
  double duration;
  start = std::clock();

  int j = 0;
  for(unsigned long i=0; i<10000000000; i++){
    j = j + 1;
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout<<"time: "<<duration<<'\n';
}

