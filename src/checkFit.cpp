#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List checkFit(arma::vec y,
              arma::mat mcmcdraws,   // MCMC iterations in rows, observations in columns.
              bool probit,
              bool doWaic,
              Rcpp::Nullable<Rcpp::NumericVector> sig = R_NilValue,
              Rcpp::Nullable<Rcpp::IntegerVector> yobs = R_NilValue){

   //*****************************************************************************
   // FUNCTION:   Calculates WAIC.
   //-----------------------------------------------------------------------------
   // INPUTS:     y           = An n-length vector of y values.  (In probit case, the latents)
   //             mcmcdraws   = A matrix of posterior draws (rows) for each obs (column). (In probit case, the latents)
   //             probit      = A boolean; if TRUE, inputs are from probit functional BART.
   //             doWaic      = A boolean; if TRUE, include WAIC calculation in output.
   //             sig         = A vector of mcmc draws for sigma. (Null in probit case)
   //             yobs        = Vector of 1/0 observed values; only populate for probit model.
   //-----------------------------------------------------------------------------
   // OUTPUT:     logdens_mcmc = Log-density (sum for all obs) for each mcmc iteration.
   //             waic         = The WAIC value for the model, if included.  Else NULL.
   //*****************************************************************************

   // Initialize structures for calculating waic.
   double p_waic;
   double lpd;
   double waic;

   // Initialize number of observations (n) and mcmc iterations (nd).
   size_t nd = mcmcdraws.n_rows;
   size_t n = mcmcdraws.n_cols;

   // Initialize vector for non-nullable sigmas, in case not using probit.
   vec sigma = zeros(nd);

   // Initialize matrix to hold individual contributions to logl for each mcmc draw
   mat ll_temp(nd,n);

   // Helper variables for probit case.
   double prob_temp;
   mat ppd(nd,n);
   vec yobs_;     // Holds non-nullable yobs, in case we use probit model.
   yobs_ = zeros(n);

   //*****************************************************************************
   //* Continuous response case.
   //*****************************************************************************

   if(probit==FALSE){

      // Convert yobs to arma vector.
      sigma = Rcpp::as<arma::vec>(sig);

      //--------------------------------------------------
      // Calculate log-density for each mcmc iteration.  (Larger is better.)
      //--------------------------------------------------

      for(size_t i=0; i<nd; i++){ // Loop through MCMC draws
         for(size_t k=0; k<n; k++){ // Loop through elements.
            ll_temp(i,k) = R::dnorm(y[k], mcmcdraws(i,k), sigma[i], TRUE);
         }
      }

      //--------------------------------------------------
      //* Calculate WAIC if indicated.
      //--------------------------------------------------

      if(doWaic){

         // Calculate p_waic and lpd, and combine to get WAIC.
         p_waic = as_scalar(sum(var(ll_temp,0,1)));
         lpd = as_scalar(sum(mean(ll_temp,1)));
         waic = -2 * (lpd - p_waic);

      } // End waic.

   } // End continuous case.

   //*****************************************************************************
   //* Probit case.
   //*****************************************************************************
   if(probit==TRUE){

      // Convert yobs to arma vector.
      yobs_ = Rcpp::as<arma::vec>(yobs);

      //--------------------------------------------------
      // Calculate log-density for each mcmc iteration.  (Larger is better.)
      //--------------------------------------------------

      for(size_t i=0; i<nd; i++){ // Loop through MCMC draws
         for(size_t k=0; k<n; k++){ // Loop through elements.
            prob_temp = R::pnorm(mcmcdraws(i,k), 0, 1, TRUE, FALSE);
            ll_temp(i,k) = R::dbinom(sum(yobs_), n, prob_temp, TRUE);
         }
      }

      //*****************************************************************************
      //* Calculate WAIC.
      //*****************************************************************************

      if(doWaic){

         // Calculate pointwise predictive density.  Loop through matrix.
         for(size_t i=0; i<nd; i++){ // Loop through MCMC draws
            for(size_t k=0; k<n; k++){ // Loop through elements.
               ppd(i,k) = R::dnorm(y[k], mcmcdraws(i,k), 1, FALSE);
            }
         }

         // Calculate p_waic and lpd, and combine to get WAIC.
         p_waic = as_scalar(sum(var(log(ppd),0,1)));
         lpd = as_scalar(sum(mean(ppd,1)));
         waic = -2 * (lpd - p_waic);

      } // End probit waic.

   } // End probit case.

   //*****************************************************************************
   // Return list of function outputs.
   //*****************************************************************************
   return(List::create(_["waic"]      = waic,      // WAIC
                       _["ll_mcmc"]   = ll_temp   // Log-loss for each mcmc iteration.
   ));
}

