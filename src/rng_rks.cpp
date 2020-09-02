#include <cmath>
#include <math.h>
#include "rng_rks.h"
#include <RcppArmadillo.h>

// Define pi.
#define M_PI        3.141592653589793238462643383280

// Generate random variate from  Kolmogorov-Smirnov distribution.
// Devroye, p. 161 - used as example of alternating series method.
// Alternating series is a clever way to do accept/reject sampling.
// See Devroy, p. 151 for the basic idea.
//
// This is slightly different than technique he uses for Jacobi paper.
//
// Here we assume density is written in either case as
// f(x) = c h(x) (1 - a_1(x) + a_2(x) + ...)
//
// Left case: (index starts at n=0)
// Devroye breaks up sum on p. 161 to get alternating sum.
// a_n = 4x^2/\pi^2 \exp(- (n^2-1) \pi^2 / (8x^2) ),   n odd  (subtract)
// a_n = (n+1)^2 \exp (- ((n+1)^2-1) \pi^2 / (8x^2) ), n even (add)
// ch  = sqrt{2\pi} \pi^2 / (4x^4) \exp(-\pi^2 / (8x^2))
//
// Right case: (index starts at n=0)
// a_n = (n+1)^2 \exp(-2x^2 ((n+1)^2-1))
// ch  = 8xe^{-2x^2}
//
// Method for generating proposals come from Lemmas 5.2 and 5.3.
//
// Just copying Devroye:

// [[Rcpp::export]]
double drawleft(double tprime) {

   bool accept = FALSE;

   // Proposal - Truncated Gamma.
   double G = 0;

   while (!accept) {
      double e1 = R::rexp(1);
      double e2 = R::rexp(1);
      e1 = e1 / (1 - 0.5 / tprime);
      e2 = 2 * e2;
      G = tprime + e1;
      accept = (e1*e1 <= tprime * e2 * (G+tprime));

      if(!accept){
         accept = (G/tprime - 1 - log(G/tprime) <= e2);
      }
   }

   // Devroye is clever.
   double X = M_PI / sqrt(8 * G);
   double W = 0.0;
   double Z = 0.5 / G;
   double P = exp(-G); // exp(- \pi^2 / (8 X^2) ).
   double n = 1;
   double Q = 1;
   double U = R::runif(0,1);

   bool go = TRUE;
   while(go) {
      W = W + Z * Q;
      if (U >= W){
         return(X);
      }
      n = n + 2;
      Q = pow(P,n*n-1);
      W = W - n*n * Q;
      go = U >= W;
   }
}

// [[Rcpp::export]]
double drawright(double trunc){
   bool my_bool = TRUE;
   while (my_bool==TRUE) {

      double E = R::rexp(1.0);
      double U = R::runif(0,1.0);
      double X = sqrt(trunc*trunc + 0.5*E);

      double W = 0;
      double n = 1;
      double Z = exp(-2 * X*X);

      bool go = TRUE;
      while (go==TRUE) {
         n = n + 1;
         W = W + n*n * pow(Z,n*n-1);
         if (U >= W){
            return(X);
         }
         n = n + 1;
         W = W - n*n * pow(Z,n*n-1);
         go = U >= W;
      }

   }

}

// [[Rcpp::export]]
double rks1() {
   double TRUNC = 0.75;
   double tprime = M_PI * M_PI / (8 * TRUNC*TRUNC);

// 0.373 = pks(0.75) - Differs from PG method.
   double p = 0.373;
   double q = 1-p;

   double X = 0;
   double u_draw = R::runif(0,1);
   if (u_draw < p/(p+q)){
      X = drawleft(tprime);
   } else{
      X = drawright(TRUNC);
   }

   return(X);
}

// [[Rcpp::export]]
Rcpp::NumericVector rks_cpp(int n){
   Rcpp::NumericVector X(n);
   for(int i=0; i<n; i++){
      X[i] = rks1();
   }
   return X;
}
