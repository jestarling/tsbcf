#include "rng_rks.h"
#include <RcppArmadillo.h>
#include "rng.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

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

//-------------------------------------------------
// Functions to support rejection sampling the
// full conditional for the lambda_i ~ KS() random variable.
//-------------------------------------------------

// This function isn't necessary; just for fun, to show how to draw GIG from cpp.
// [[Rcpp::export]]
Rcpp::NumericVector rinvgauss_cpp(size_t n, double mean, double shape, double dispersion){

   // calling rnorm()
   Function f("rinvgauss");

   // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
   return f(n, Named("mean")=mean, _["shape"]=shape, _["dispersion"]=dispersion);
}

// These three functions are from Holmes & Held (2006), Appendix A4.
// https://projecteuclid.org/download/pdf_1/euclid.ba/1340371078

// [[Rcpp::export]]
double rks_posterior(double r2){
   // Note: r^2 = (z_i -fit_i)^2
   double r = sqrt(r2);

   bool ok = false;

   while(ok==false){

      // To begin, we must draw a sample from the rejection sampling density.
      // Using w instead of y to avoid confusion with y in main tsbcfLogic.cpp file.
      double w = randn(); // w ~ N(0,1)
      w = w*w;
      w = 1 + (w - sqrt(w * (4*r+w))) / (2*r);
      double u = randu();

      double lambda;
      if(u <= (1/(1+w))){
         lambda = r/w;
      } else{
         lambda = r*w;
      }

      if(lambda<0){
         Rcpp::Rcout << "Lambda is negative: " << lambda << endl;
         Rcpp::Rcout << "r: " << r << ", w: " << w << endl;
      }

      // Now, the rejection sampling part.
      // Now lambda ~ GIG(.5, 1, r^2)
      u = randu();

      if(lambda>(4/3)){
         ok = rightmost_interval(u,lambda);
         return(lambda);
      } else{
         ok = leftmost_interval(u,lambda);
         return(lambda);
      }
   }
}

// [[Rcpp::export]]
double rightmost_interval(double u, double lambda){

   double z = 1;
   double x = exp(-.5*lambda);

   // Repeat
   double j = 0;

   bool time_to_stop = false;
   bool ok = false;

   while(!time_to_stop) {

      // Squeezing.
      j = j+1;
      z = z - (j+1) * (j+1) * pow(x,(j+1)*(j+1)-1);

      if(z>0){
         ok = true;
         time_to_stop=true;
         return(ok);
      }

      j = j+1;
      z = z - (j+1) * (j+1) * pow(x,(j+1)*(j+1)-1);

      if(z<u){
         ok = false;
         time_to_stop=true;
         return(ok);
      }
   }
}

// [[Rcpp::export]]
double leftmost_interval(double u, double lambda){

   double h = .5*log(2) + 2.5*(log(M_PI)) - 2.5*log(lambda) - M_PI*M_PI/(2*lambda) + .5*lambda;
   double lu = log(u);
   double z = 1;
   double x = exp(-M_PI*M_PI/(2*lambda));
   double k = lambda/(M_PI*M_PI);
   double j = 0;

   bool time_to_stop = false;
   bool ok = false;

   while(!time_to_stop) {

      // Squeezing.
      j = j+1;
      z = z - k*pow(x,j*j-1);

      if(h + log(z) > lu){
         ok = true;
         time_to_stop = true;
         return(ok);
      }

      j = j+1;
      z = z + (j+1)*(j+1)*pow(x,(j+1)*(j+1)-1);

      if(h+log(z)<lu){
         ok = false;
         time_to_stop = true;
         return(ok);
      }
   }
}




