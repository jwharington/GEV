
#include <iostream>
#include <vector>
extern "C" {
#include "as215.h"
}

int main() {

  std::vector<double> x;
  double para[3] = {0,0,0};

  // https://www.jstor.org/stable/2347483?seq=1

  int n=0;
  while (1) {
    double _x;
    std::cin >> _x;
    if (std::cin.eof()) {
      break;
    }
    x.push_back(_x);
    n++;
  }
  printf("read %d\n", n);

  while (n>2) {
    int monit = 0; // 1 for feedback
    int ifault;
    double vcov[6];
    mlegev_(x.data(),&n,para,vcov,&monit,&ifault);
    if (ifault==0) {
      printf("parameters %g %g %g with %d\n", para[0], para[1], para[2], n);
      break;
    }
    n-= n/100;
  }
  // returns: mu, sigma, xi
}

