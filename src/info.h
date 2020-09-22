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
   arma::vec sigma_het; // For logit formulation, with heterogeneous sigma_i.

   //----------------------------
   // Original pinfo constructor; no time-series length argument.
   //----------------------------
   pinfo() {pbd=1.0; pb=.5; alpha=.95; beta=.5;  sigma=1.0;}

   //----------------------------
   // New properties for tsbcf.
   //----------------------------

   arma::vec mu0;       // Vector of prior means; length tlen.
   double ls;           // length-scale parameter for the Squared Exponential kernel.
   double tau;          // scalar sd for the Squared Exponential kernel. (Does not matrix part of SE kernel.)
   mat K;               // Carries unscaled part of SE kernel.
   mat Prec;            // Precision matrix for mu_l ~iid N(mu0, inv(Prec))
   mat Prec0;           // Prior precision matrix.

   // Note: SE kernel = (tau^2/2) * exp(.5 (( t-t')/l)^2)

   //----------------------------
   // tsbcf pinfo constructor
   //----------------------------

   pinfo(size_t tlen) {
      pbd=1.0;
      pb=.5;
      alpha=.95;
      beta=0.5;
      sigma=1.0;
      mu0 = zeros(tlen);
      ls=1.0;
      tau=1.0;
      Prec = eye(tlen,tlen);
      }

   // Constructor for logit - provides n, to initialize sigma_het vector.
   pinfo(size_t tlen, size_t n) {
      pbd=1.0;
      pb=.5;
      alpha=.95;
      beta=0.5;
      sigma_het=ones(n);
      mu0 = zeros(tlen);
      ls=1.0;
      tau=1.0;
      Prec = eye(tlen,tlen);
   }

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
