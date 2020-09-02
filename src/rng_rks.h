#ifndef RNG_H
#define RNG_H
#include <RcppArmadillo.h>

using std::vector;


double drawleft(double tprime);
double drawright(double trunc);
double rks1();
double rks(int n=1);
#endif // RNG_H

