#include <RcppArmadillo.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector rpgmix(int n, double a, int m) {
   NumericVector out(n); out.fill(0.0);
   double loga = log(a);
   double t;
   for(int i=0; i<n; ++i) {
     for(int j=0; j<m; ++j) {
       t = log(R::rgamma(a, 1.0)) - loga;
       if(R::runif(0.0, 1.0)<0.5) t = -t;
       out(i) += t;
     }
   }
   return out;
}
