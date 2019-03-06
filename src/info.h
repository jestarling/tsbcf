#ifndef GUARD_info_h
#define GUARD_info_h

#include <RcppArmadillo.h>
using namespace arma;

//============================================================
//data
//============================================================

class dinfo {
public:

   // Methods
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j).  *x is first element of the data - is a pointer.
   double *y; // ith y is *(y+i) or y[i]

   double *t; // ith t is *(t+i) or t[i].  (*t)[i] to reference value.  (Points to first element of t vector.)
   size_t tlen; // size_t T;  // Length of vector t_ref; total number of time points.
   vec tref; // Vector of unique t values.

   size_t nt; // Number of treated

   // Constructor
   dinfo() {p=0;n=0;x=0;y=0; tlen = 1; t=0;}

};

//============================================================
//prior and mcmc
//============================================================

class pinfo
{
public:

   //----------------------------
   // Declare properties.
   //----------------------------

   //mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
   //prior info
   double alpha;
   double beta;
   //sigma
   double sigma;

   //----------------------------
   // Original pinfo constructor; no time-series length argument.
   //----------------------------
   pinfo() {pbd=1.0; pb=.5; alpha=.95; beta=.5;  sigma=1.0;}

   //----------------------------
   // new properties for functional bart.
   //----------------------------

   arma::vec mu0;       // Vector of prior means; length tlen.
   double ls;           // length-scale parameter for the Squared Exponential kernel.
   double var;          // variance parameter for the Squared Exponential kernel.
   mat Sigma0;          // Prior cov matrix for mu_l ~iid N(mu0, Sigma0 = K^-1)
   mat Prec0;           // Prior precision matrix for mu_l, ie inv(Sigma0).

   // For augmented model to induce C+(0,sig2) hyperprior on var.
   // y = eta * f(x) + e, e ~ N(0,sig2)
   // eta | gamma ~ N(0,gamma2)
   double eta;
   double gamma;

   // //----------------------------
   // // new properties for causal functional bart.
   // //----------------------------
   // double eta_m;
   // double eta_a;
   // double gamma_m;
   // double gamma_a;
   // double beta_trt;   // Holds separate beta to regularize the alpha(x,pix)*zi trees.

   //----------------------------
   // Functional bart pinfo constructor
   //----------------------------

   pinfo(size_t tlen) { pbd=1.0; pb=.5; alpha=.95; beta=0.5; sigma=1.0;
   mu0 = zeros(tlen);
   ls=1.0; var=.005;  // Default for var is 1/m where m = 200 by default
   eta = 1.0; gamma = 1.0;
   Sigma0 = eye(tlen,tlen);
   Prec0 = eye(tlen,tlen);}

   // //----------------------------
   // // Causal functional bart pinfo constructor
   // //----------------------------
   //
   // pinfo(size_t tlen, double beta_t) { pbd=1.0; pb=.5; alpha=.95; beta=2; sigma=1.0;
   // beta_trt = beta_t;
   // mu0 = zeros(tlen);
   // ls=1.0; var=.005;  // Default for var is 1/m where m = 200 by default
   // eta_m = 1.0; eta_a = 1.0;
   // gamma_m = 1.0; gamma_a = 1.0;
   // Sigma0 = eye(tlen,tlen);}

};

//============================================================
//sufficient statistics for 1 node
//============================================================

class sinfo
{
public:

   double n0; //unweighted sample counts for het case.
   double n;
   double sy;
   double sy2;
   vec n0_vec; // unweighted sample counts for het case at each t.
   vec n_vec;
   vec sy_vec;

   // Original constructor; no time series length argument.
   sinfo() {n0=0.0;n=0;sy=0.0;sy2=0.0; n0_vec = zeros(1); n_vec = zeros(1); sy_vec = zeros(1);}

   // New constructor, with time series length argument.
   sinfo(size_t tlen) {n0=0.0;n=0;sy=0.0;sy2=0.0; n0_vec = zeros(tlen); n_vec = zeros(tlen); sy_vec = zeros(tlen);}

};

#endif
