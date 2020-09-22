#ifndef RNG_H_KS
#define RNG_H_KS
#include <RcppArmadillo.h>

using std::vector;

// Define pi.
#define M_PI        3.141592653589793238462643383280


double drawleft(double tprime);
double drawright(double trunc);
double rks1();
Rcpp::NumericVector rinvgauss_cpp(size_t n, double mean, double shape, double dispersion);
double rks_posterior(double r);
double rightmost_interval(double u, double lambda);
double leftmost_interval(double u, double lambda);

#endif // RNG_H

