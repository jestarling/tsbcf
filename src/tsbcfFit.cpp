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
List tsbcfFit(arma::vec y,
              arma::vec z,
              arma::vec zpred,
              arma::vec tgt,
              arma::vec tpred,
              arma::vec x_con,
              arma::vec x_mod,
              arma::vec xpred_con,
              arma::vec xpred_mod,
              List xinfo_list_con,
              List xinfo_list_mod,
              arma::vec trt_init,
              int nburn, int nsim, int ntree_con=200, int ntree_mod=50,
              double lambda=-999, double sigq=.9, double sighat=-999, double nu=3,
              double base_con=.95, double power_con=2.0,
              double base_mod=.25, double power_mod=3.0,
              double ecross_con=1, double ecross_mod=1,
              double con_sd=1, double mod_sd=1,
              bool use_muscale = false,
              bool use_tauscale = false,
              CharacterVector treef_name_="tsbcf_trees.txt",
              bool save_trees = false,
              bool silent_mode = false)
{

   // Begin output stream.
   std::string treef_name = as<std::string>(treef_name_);
   std::ofstream treef(treef_name.c_str());

   // If save_trees = false, set badbit so no output generated.
   if(save_trees==false){
      treef.setstate(std::ios_base::badbit);
   }

   // Random number generator, used in all draws.
   RNGScope scope;
   RNG gen;

   if(silent_mode==false){ Rcout << "\n*****Into bart main\n"; }

   //*****************************************************************************
   // Sample sizes for trt obs and full data set.
   //*****************************************************************************

   // In-sample n and ntrt.
   size_t n = y.size();
   size_t ntrt = sum(z);

   // Out-of-sample n and ntrt.
   size_t n_pred = tpred.size();
   size_t ntrt_pred = sum(zpred);

   //*****************************************************************************
   // Read, format info about times
   //*****************************************************************************

   vec tref = unique(join_cols(tgt,tpred));  // Vector of unique time points, sorted in asc. order.
   int tlen = tref.n_elem;                   // Number of unique time points (T).

   if(silent_mode==false){
      Rcout << "unique times: " << endl << tref << endl;
      Rcout << "Number of total time points T: " << tlen << endl;
   }

   //------------------------------------------------------
   // Calculate alpha_t (ybar_t) vector, for conditional-on-time mean-centering data.
   //------------------------------------------------------

   vec alpha_t = zeros(tlen);    // Holds ybar at each t.
   vec n_t = zeros(tlen);        // Holds sample sizes at each t.
   uvec idx;                     // Holds indices where ti = t, for calculating ybar_t.

   // Iterate through time points.
   for(int i = 0; i < tlen; i++){

      // All observations.
      idx = find(tgt == tref(i));
      n_t(i) = idx.size();

      if(n_t(i) > 0){
         alpha_t(i) = sum(y(idx)) / n_t(i);
      }
   }

   // Create n-length vector with corresponding ybar_i for each y_i.
   // Will let us subtract off time-means and add back later.
   vec alpha_ti(y.size());
   idx.reset();

   for(size_t i = 0; i < y.size(); i++){
      idx = find(tref == tgt(i));
      alpha_ti(i) = as_scalar(alpha_t(idx));
   }

   //*****************************************************************************
   // Format y, and calculate sufficient statistics for y, for initializing
   //*****************************************************************************

   //------------------------------------------------------
   // Mean-center data, conditional on the mean for each time point.
   //------------------------------------------------------
   y = y - alpha_ti;

   //------------------------------------------------------------------------
   // Scalar versions of sufficient stats for all y.
   //------------------------------------------------------------------------

   sinfo allys;
   allys.n = n;
   allys.sy = sum(y);
   allys.sy2 = sum(y % y);

   //------------------------------------------------------------------------
   // Vector sufficient stats for y_t's (T-length time vector).
   //------------------------------------------------------------------------

   vec ybar_vec = zeros(tlen);
   allys.n_vec = n_t;
   allys.sy_vec = zeros(tlen);

   // Iterate through time points; calculate ybar_t.
   for(int k = 0;k<tlen;k++){
      if(n_t(k)>0){
         idx = find(tref == tgt(k));
         allys.sy_vec(k) = sum(y(idx));
         ybar_vec[k] = allys.sy_vec(k) / n_t[k];
      }
   }


   //*****************************************************************************
   //* Read, format X's: x_con and x_mod for in-samp, and xpred_con, xpred_mod for oos.
   //* The n*p numbers for x are stored as the p for first obs, then p for second, and so on.
   //*****************************************************************************

   size_t p_con = x_con.size()/n;   // Number of prognostic covariates
   size_t p_mod = x_mod.size()/n;   // Number of treatment covariates.

   size_t p_con_pred = xpred_con.size()/n_pred; // Number of oos prognostic covariates.
   size_t p_mod_pred = xpred_mod.size()/n_pred; // Number of oos treatment covariates.

   // Check correct number of covariates.
   double npred_con = xpred_con.size()/p_con;   // Number of covariates for prognostic.
   double npred_mod = xpred_mod.size()/p_mod;   // Number of covariates for treatment.

   if(xpred_con.size() != npred_con*p_con) Rcout << "error, wrong number of elements in prediction data set\n";
   if(xpred_mod.size() != npred_mod*p_mod) Rcout << "error, wrong number of elements in prediction data set for treatment effects\n";

   //------------------------------------------------------------------------
   //x cutpoints
   //------------------------------------------------------------------------

   xinfo xi_con;
   xinfo xi_mod;

   // prognostic
   xi_con.resize(p_con);
   for(int i=0; i<p_con; ++i) {
      NumericVector tmp = xinfo_list_con[i];
      std::vector<double> tmp2;
      for(size_t j=0; j<tmp.size(); ++j) {
         tmp2.push_back(tmp[j]);
      }
      xi_con[i] = tmp2;
   }

   // treatment
   xi_mod.resize(p_mod);
   for(int i=0; i<p_mod; ++i) {
      NumericVector tmp = xinfo_list_mod[i];
      std::vector<double> tmp2;
      for(size_t j=0; j<tmp.size(); ++j) {
         tmp2.push_back(tmp[j]);
      }
      xi_mod[i] = tmp2;
   }


   //*****************************************************************************
   // Setup for MCMC
   //*****************************************************************************

   //------------------------------------------------------------------------
   // PX scale parameter for tau: tau(x,t) = bscale * tau_0(x,t) where
   // b1 ~ N(mod_sd/2, 1/2) and b0 ~ N(-mod_sd/2, 1/2) s.t. (b1-b0) ~ N(mod_sd, 1)
   //------------------------------------------------------------------------

   double bscale1 = 0.5;   // Initialize value for b1.
   double bscale0 = -0.5;  // Initialize value for b0.
   double bscale1_mean = mod_sd/2;
   double bscale0_mean = -mod_sd/2;
   double bscale_prec = 2;

   //------------------------------------------------------------------------
   // Initialize trees.
   //------------------------------------------------------------------------

   // mu(x,t,pi) for prognostic trees. if you sum the fit over the trees you get the fit.
   std::vector<tree> t_con(ntree_con);
   for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar_vec / ntree_con);

   // tau(x,t) for treatment trees
   std::vector<tree> t_mod(ntree_mod);
   for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init / ntree_mod);

   //------------------------------------------------------------------------
   // Prior for prognostic trees.
   //------------------------------------------------------------------------
   pinfo pi_con(tlen);

   pi_con.pbd = 1.0;          //prob of birth/death move
   pi_con.pb = .5;            //prob of birth given  birth/death
   pi_con.alpha = base_con;   //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi_con.beta = power_con;   //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node

   // Sigma.  Use user-provided sighat, unless null (-999); then use sd(y).  Gets adjusted later for backfitting: sig^2 / eta_con^2.
   if(sighat!=-999){
      pi_con.sigma = sighat;
   } else{
      pi_con.sigma = stddev(y);
   }

   // Use override for lambda if provided by user.  Else calculate using sighat and sigq.
   double qchi = 0;
   if(lambda==-999){
      qchi = R::qchisq(1-sigq, nu, true, false);
      lambda = (sighat * sighat * qchi) / nu;
   }

   // Initialize new time-related prior components.
   pi_con.mu0 = zeros(tlen);
   pi_con.ls = tlen / (PI * ecross_con);

   // If use_muscale=false, then var(f(x,t)) = con_sd^2 / m is variance of Cov matrix.
   // If use_muscale=true, then pi.var = 1/m, and C+ prior is induced with median con_sd via eta.
   pi_con.var = 1 / static_cast<double>(ntree_con);
   if(use_muscale==false) pi_con.var = pi_con.var * con_sd * con_sd;

   pi_con.Sigma0 = cov_se(tref, tref, pi_con.ls, pi_con.var);
   pi_con.Prec0 = pi_con.Sigma0.i(); // Prior precision matrix.
   pi_con.eta = 1;
   pi_con.gamma = 1;

   //------------------------------------------------------------------------
   // Prior for treatment trees.
   //------------------------------------------------------------------------
   pinfo pi_mod(tlen);

   pi_mod.pbd = 1.0;                // prob of birth/death move
   pi_mod.pb = .5;                  // prob of birth given  birth/death
   pi_mod.alpha = base_mod;         // prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi_mod.beta = power_mod;         // 3 for treatment trees.
   pi_mod.sigma = pi_con.sigma;     // Same initial value as pi_con.sigma.  Gets adjusted later for backfitting: sig^2 / eta_mod^2.

   // Initialize new time-related prior components.
   pi_mod.mu0 = zeros(tlen);
   pi_mod.ls = tlen / (PI * ecross_mod);

   // If use_tauscale=false, then var(f(x,t)) = con_sd^2 / m is variance of Cov matrix.
   // If use_tauscale=true, then pi.var = 1/m, and C+ prior is induced with median con_sd via eta.
   pi_mod.var = 1 / static_cast<double>(ntree_con);
   if(use_tauscale==false) pi_mod.var = pi_mod.var * mod_sd * mod_sd;

   pi_mod.Sigma0 = cov_se(tref, tref, pi_mod.ls, pi_mod.var);
   pi_mod.Prec0 = pi_mod.Sigma0.i(); // Prior precision matrix.


   //------------------------------------------------------------------------
   // dinfo for control function mu(x,t)
   //------------------------------------------------------------------------
   double* allfit_con = new double[n];    // Holds sum of fit of all trees for each observation.
   double* r_con = new double[n];         // Holds residuals.  y - (allfit - ftemp) = y - allfit + ftemp.

   // Initialize allfit_con to the ybar_t for each obs' corresponding time t value.
   for (int i=0; i<n; ++i){
      idx = find(tref==tgt(i));
      allfit_con[i] = as_scalar(ybar_vec(idx));
   }

   dinfo di_con;
   di_con.n=n; di_con.p=p_con; di_con.x = &x_con[0]; di_con.y=r_con; //the "y" for each draw will be the residual
   di_con.tlen = tlen; di_con.tref = tref; di_con.t = &tgt[0];

   //------------------------------------------------------------------------
   // dinfo for treatment function tau(x,t)
   //------------------------------------------------------------------------
   double* allfit_mod = new double[n];    // Holds sum of fit of all trees for each observation.
   double* r_mod = new double[n];         // Holds residuals.  y - (allfit - ftemp) = y - allfit + ftemp.

   // Initialize allfit to the trt_init_t = (ybar_trt_t - ybar_ctrl_t) for each obs' corresponding time t value.
   for (int i=0; i<n; ++i){
      idx = find(tref==tgt(i));
      allfit_mod[i] = as_scalar( (z[i]*bscale1 + (1-z[i])*bscale0) * trt_init(idx));
   }

   dinfo di_mod;
   di_mod.n=n; di_mod.p=p_mod; di_mod.x = &x_mod[0]; di_mod.y=r_mod; //the "y" for each draw will be the residual
   di_mod.tlen = tlen; di_mod.tref = tref; di_mod.t = &tgt[0];

   //------------------------------------------------------------------------
   // dinfo for prognostic function mu(x,t) out-of-sample
   //------------------------------------------------------------------------
   dinfo di_con_pred;
   di_con_pred.n = 0;

   if(n_pred>0){
      di_con_pred.n = n_pred; di_con_pred.p = p_con_pred; di_con_pred.x = &xpred_con[0]; di_con_pred.y=0; //there are no y's.
      di_con_pred.tref = unique(tpred); di_con_pred.tlen = tlen; di_con_pred.t = &tpred[0];
   }

   //------------------------------------------------------------------------
   // dinfo for treatment function tau(x,t) out-of-sample
   //------------------------------------------------------------------------
   dinfo di_mod_pred;
   di_mod_pred.n = 0;

   if(n_pred>0){
      di_mod_pred.n = n_pred;
      di_mod_pred.p = p_mod_pred;
      di_mod_pred.x = &xpred_mod[0];
      di_mod_pred.y=0; //there are no y's.
      di_mod_pred.tref = unique(tpred);
      di_mod_pred.tlen = tlen;
      di_mod_pred.t = &tpred[0];
   }

   // TROUBLESHOOTING:
   //  Rcout << "Completed di section" << endl;
   //  Rcout << "Checking pointer: di.t" << endl << di.t;
   //  Rcout << "Checking pointer: di.t value" << endl << di.t[0];

   //------------------------------------------------------------------------
   // storage for entire fit (yhat)
   //------------------------------------------------------------------------
   double* allfit = new double[n]; //yhat

   for(size_t k=0;k<n;k++){
      allfit[k] = allfit_con[k] + allfit_mod[k];
   }

   double* ftemp = new double[n];   // Holds fit of current tree.
   double* precs = new double[n];   // Temp storage for conditional ''precisions'' in heteroskedastic updates.

   //------------------------------------------------------------------------
   //storage for full conditionals and associated quantities.
   //------------------------------------------------------------------------

   // For sigma draw
   NumericVector sigma_post(nsim);

   // For prognostic scale.
   NumericVector eta_post_con(nsim);
   NumericVector gamma_post_con(nsim); // SD for prognostic prior's scale.

   // For components of treatment scale.
   NumericVector bscale1_post(nsim);
   NumericVector bscale0_post(nsim);

   // For holding posterior draws for in-sample and out-of-sample data.
   NumericMatrix yhat_post(nsim,n);
   NumericMatrix con_post(nsim,n);
   NumericMatrix mod_post(nsim,n);
   NumericMatrix yhat_est_post(nsim,n_pred);
   NumericMatrix con_est_post(nsim,n_pred);
   NumericMatrix mod_est_post(nsim,n_pred);

   // For holding MH alphas for prognostic and treatment trees.
   NumericMatrix alpha_con (nsim+nburn,ntree_con);
   NumericMatrix alpha_mod (nsim+nburn,ntree_mod);

   //------------------------------------------------------------------------
   //save stuff to tree file
   //------------------------------------------------------------------------
   double thin = 1;

   treef << xi_con << endl;         //cutpoints
   treef << xi_mod << endl;         //cutpoints
   treef << ntree_con << endl;      //number of trees
   treef << ntree_mod << endl;      //number of trees
   treef << p_con << endl;          //dimension of x's
   treef << p_mod << endl;          //dimension of x's
   treef << (int)(nsim/thin) << endl;

   //*****************************************************************************
   //* MCMC
   //*****************************************************************************
   if(silent_mode==false){ Rcout << "\nMCMC:\n"; }
   time_t tp;
   int time1 = time(&tp);

   for(size_t i=0;i<(nsim+nburn);i++) {

      // Progress report.
      if(silent_mode==false){ if(i%50==0) cout << "MCMC ITERATION: " << i << " of " << nsim+nburn <<  endl; }

      //------------------------------------------------------------------------
      //draw prognostic trees for [alpha_mu_t + eta_mu * mu(x,t,pihat)]
      //------------------------------------------------------------------------

      // Loop through ntree_con prognostic trees.
      for(size_t j=0;j<ntree_con;j++) {

         // Grab fit of current tree.
         fit(t_con[j], xi_con, di_con, ftemp);

         // Loop through observations.
         for(size_t k=0;k<n;k++) {

            // Check for NA in tree fits.
            if(ftemp[k] != ftemp[k]){
               Rcout << t_con[j] << endl;
               stop("nan in ftemp");
            }

            // Subtract out jth tree from fits and calculate residual.
            allfit[k] = allfit[k] - pi_con.eta*ftemp[k];
            allfit_con[k] = allfit_con[k] - pi_con.eta*ftemp[k];
            r_con[k] = (y[k]-allfit[k])/pi_con.eta;
            //r_con[k] = (y[k]-allfit_con[k]-allfit_mod[k])/pi_con.eta;

            // Check for NA in residual.
            if(r_con[k] != r_con[k]){
               stop("nan in con residual");
            }
         } // End loop over n observations for jth tree.

         // Birth-death proposal and drawing new leaf mu's.
         alpha_con(i,j) = bd(t_con[j],xi_con,di_con,pi_con,gen);
         drmu(t_con[j],xi_con,di_con,pi_con,gen);

         // Update the current tree, and add back the updated jth trees subtracted from each sum.
         fit(t_con[j],xi_con,di_con,ftemp);

         for(size_t k=0;k<n;k++){
            allfit[k] += pi_con.eta*ftemp[k];
            allfit_con[k] += pi_con.eta*ftemp[k];
         }

      } // End prognostic tree loop.


      //------------------------------------------------------------------------
      //draw treatment trees for [(alpha_tau_t + eta_tau * mu(x,t,pihat) ) * bscale]
      //------------------------------------------------------------------------

      // Set up temporary 'precisions' for heteroskedastic updates.
      for(size_t k=0;k<ntrt;k++){
         precs[k] = (bscale1*bscale1)/(sighat*sighat);
      }
      for(size_t k=ntrt;k<n;k++){
         precs[k] = (bscale0*bscale0)/(sighat*sighat);
      }

      // Loop through ntree_mod treatment trees.
      for(size_t j=0;j<ntree_mod;j++) {

         // Grab fit of current tree.
         fit(t_mod[j],xi_mod,di_mod,ftemp);

         // Loop through observations.
         for(size_t k=0;k<n;k++) {

            // Check for NA in tree fits.
            if(ftemp[k] != ftemp[k]) {
               Rcout << t_mod[j] << endl;
               stop("nan in ftemp");
            }

            // Subtract out jth tree from fits and calculate residual.
            double bscale = (k<ntrt) ? bscale1 : bscale0;
            allfit[k] -= bscale*ftemp[k];
            allfit_mod[k] -= bscale*ftemp[k];
            r_mod[k] = (y[k]-allfit[k])/bscale;

            // Check for NA in residual.
            if(r_mod[k] != r_mod[k]){
               stop("nan in mod residual");
            }

         } // End loop over n observations for jth tree.

         // Birth-death proposal and drawing new leaf mu's.
         alpha_mod(i,j) = bdhet(t_mod[j],xi_mod,di_mod,precs,pi_mod,gen);
         drmuhet(t_mod[j],xi_mod,di_mod,precs,pi_mod,gen);

         // Update the current tree, and add back the updated jth trees subtracted from each sum (trt, then control, for scaling).
         fit(t_mod[j],xi_mod,di_mod,ftemp);

         for(size_t k=0;k<ntrt;k++){
            allfit[k] += bscale1*ftemp[k];
            allfit_mod[k] += bscale1*ftemp[k];
         }
         for(size_t k=ntrt;k<n;k++){
            allfit[k] += bscale0*ftemp[k];
            allfit_mod[k] += bscale0*ftemp[k];
         }

      } // End treatment tree loop.

      //------------------------------------------------------------------------
      // Draw eta_con (scale for prognostic tree fits).
      //------------------------------------------------------------------------

      if(use_muscale==true){

         double ww = 0.0;  // Tracks sum of w*w
         double rw = 0.0;  // Tracks sum of r*w

         for(size_t k=0;k<n;k++) {
            double r = (y[k] - allfit_mod[k]) * pi_con.eta / allfit_con[k]; // Holds r_i
            double w = allfit_con[k] * allfit_con[k] / (pi_con.eta * pi_con.eta); // Holds tau(x,t)^2_i

            ww += w*w;
            rw += r*w*w;
         }

         // Calculate variance and mean for eta full conditional. eta ~ N(con_sd, gamma2) to make con_sd the C+ prior's median.
         double eta_var = 1 / (1/(pi_con.gamma*pi_con.gamma) + ww/(sighat*sighat));
         double eta_mean = eta_var * (con_sd/(pi_con.gamma*pi_con.gamma) + rw/(sighat*sighat));

         // Draw normal full conditional for eta.
         double eta_old = pi_con.eta; // Save previous eta before drawing new one, for adjusting scaling.
         pi_con.eta = eta_mean + gen.normal(0.,1.) * sqrt(eta_var);

         // Update prognostic fits to have new pi.eta scaling.
         for(size_t k=0; k<n; ++k) {
            allfit_con[k] = allfit_con[k] * pi_con.eta / eta_old;
         }

      } else{
         pi_con.eta = 1;
      }

      //------------------------------------------------------------------------
      // draw b1 and b0.
      //------------------------------------------------------------------------

      if(use_tauscale==true){

         // Calculate sums required for means and variances.
         double ww0 = 0.0, ww1 = 0.0; // holds precisions.
         double rw0 = 0.0, rw1 = 0.0; // holds resids*precisions.
         double s2 = sighat * sighat;

         for(size_t k=0;k<n;k++){
            double bscale = (k<ntrt) ? bscale1 : bscale0;
            double w = (allfit_mod[k] * allfit_mod[k])/ (bscale * bscale);

            if(w!=w){
               Rcout << "w: " << w << endl;
               stop("");
            }

            double r = (y[k] - allfit_con[k]) * bscale / allfit_mod[k];

            if(r!=r){
               Rcout << "bscale " << k << " r " << r << endl;
               stop("");
            }

            if(k<ntrt){
               ww1 += w;
               rw1 += r*w;
            } else{
               ww0 += w;
               rw0 += r*w;
            }
         }

         // Calculate mean and variance for b0 and b1 full conditionals, and draw.
         double bscale1_old = bscale1;
         double bscale_var = 1 / (ww1/s2 + bscale_prec);
         double bscale_mean = bscale_var * (rw1/s2 + bscale_prec*bscale1_mean);
         bscale1 = gen.normal(0.,1.) * sqrt(bscale_var) + bscale_mean;

         double bscale0_old = bscale0;
         bscale_var = 1 / (ww0/s2 + bscale_prec);
         bscale_mean = bscale_var * (rw0/s2 + bscale_prec*bscale0_mean);
         bscale0 = gen.normal(0.,1.) * sqrt(bscale_var) + bscale_mean;

         // Rescale fits.
         for(size_t k=0;k<ntrt;k++){
            allfit_mod[k] = allfit_mod[k]*bscale1/bscale1_old;
         }
         for(size_t k=ntrt;k<n;k++){
            allfit_mod[k] = allfit_mod[k]*bscale0/bscale0_old;
         }

      } else{
         bscale1 = .5;
         bscale0 = -.5;
      }

      //------------------------------------------------------------------------
      // draw gamma_con (prognostic)
      //------------------------------------------------------------------------
      pi_con.gamma = sqrt((1 + ((pi_con.eta-con_sd) * (pi_con.eta-con_sd))) / gen.chi_square(2));

      //------------------------------------------------------------------------
      // Sync yhat=allfit[k] after scale updates.
      //------------------------------------------------------------------------

      // allfit_con and allfit_mod are just the trees, without scales/centers; allfit holds entire yhat.
      for(size_t k=0;k<n;k++){
         allfit[k] = allfit_con[k] + allfit_mod[k];
      }

      //------------------------------------------------------------------------
      // Draw sigma.
      //------------------------------------------------------------------------

      // Calculate residuals.
      double rss = 0.0;
      for(size_t k=0;k<n;k++) {
         double r = y[k] - allfit[k];
         rss += r*r;
      }

      // Draw sigma.
      sighat = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));

      // Update for backfitting.
      pi_con.sigma = sighat / fabs(pi_con.eta);
      pi_mod.sigma = sighat;

      //------------------------------------------------------------------------
      // Save MCMC output.
      //------------------------------------------------------------------------

      // First, need a vector for the alpha_ti values at each tpred.
      // For out of sample observations, pull in appropriate ybar_t
      // from each time point.  (Uses ybar_t from in-sample data for each time point).
      double np = zpred.size();
      vec alpha_ti_pred(np);
      for(size_t i = 0; i < np; i++){
         idx = find(tref == tpred(i));
         alpha_ti_pred(i) = as_scalar(alpha_t(idx));
      }

      // Save only if exceeded number of iterations to nburn.
      if(i>=nburn) {

         // Save the ntree_mod tree outputs.
         for(size_t j=0;j<ntree_mod;j++) treef << t_mod[j] << endl;

         // Save vectors of sigma, eta and gamma draws.
         sigma_post(i-nburn) = sighat;
         eta_post_con(i-nburn) = pi_con.eta;
         gamma_post_con(i-nburn) = pi_con.gamma;
         bscale1_post(i-nburn) = bscale1;
         bscale0_post(i-nburn) = bscale0;

         // Save in-sample treatment effects, prognostic effects, and yhat fits.
         for(size_t k=0;k<n;k++) {

            yhat_post(i-nburn,k) = allfit[k] + alpha_ti[k];
            double bscale = (k<ntrt) ? bscale1 : bscale0;
            mod_post(i-nburn,k) = (bscale1-bscale0) * allfit_mod[k]/bscale; // allfit_mod[k] = [b1*z_i + b0*(1-z_i)] * tau(x,t) == bscale * tau(x,t) --> need (b1 - b0) * tau(x,t), so divide out bscale.
            con_post(i-nburn,k) = yhat_post(i-nburn,k) - z[k]*mod_post(i-nburn,k);
         }

         // Draw out of sample treatment effects.
         if(di_mod_pred.n>0){
            for(size_t k=0;k<n_pred;k++) {
               double tau_temp = fit_i(k,t_mod,xi_mod,di_mod_pred);
               double bscale = (k<ntrt_pred) ? bscale1 : bscale0;

               mod_est_post(i-nburn,k) = (bscale1 - bscale0) * tau_temp;
               yhat_est_post(i-nburn,k) = alpha_ti_pred[k] + pi_con.eta * fit_i(k,t_con,xi_con,di_con_pred) + tau_temp*bscale;
               con_est_post(i-nburn,k) =  yhat_est_post(i-nburn,k) - zpred[k]*mod_est_post(i-nburn,k);
            }
         }

      } // End saving MCMC outputs.
   } // End MCMC loop.

   //-------------------------------------------------
   // MCMC loop time keeping
   //-------------------------------------------------

   if(silent_mode==false){   Rcout << endl << "MCMC LOOP FINISHED SUCCESSFULLY" << endl; }
   int time2 = time(&tp);
   if(silent_mode==false){ Rcout << "time for loop: " << time2 - time1 << endl; }


   //*****************************************************************************
   //* Clean up and return results.
   //*****************************************************************************

   t_con.clear();
   t_mod.clear();
   delete[] allfit_con;
   delete[] allfit_mod;
   delete[] r_con;
   delete[] r_mod;
   delete[] ftemp;
   treef.close();

   // If save_tree = false, delete tree file.
   if(save_trees==false){
      std::remove(treef_name.c_str());
   }

   return(List::create(_["yhat"]       = yhat_post,         // (nsim x n) matrix of in-sample mcmc yhats.
                       _["yhat_oos"]   = yhat_est_post,     // (nsim x npred) matrix of out-of-sample yhats.
                       _["mu"]         = con_post,          // (nsim x n) matrix of in-sample mcmc prognostic effects.
                       _["mu_oos"]     = con_est_post,      // (nsim x npred) matrix of out-of-sample mcmc prognostic effects.
                       _["tau"]        = mod_post,          // (nsim x n) matrix of in-sample mcmc treatment effects.
                       _["tau_oos"]    = mod_est_post,      // (nsim x npred) matrix of out-of-sample mcmc treatment effects.
                       _["sigma"]      = sigma_post,        // nsim-length vector of sigma draws.
                       _["bscale1"]    = bscale1_post,      // nsim-length vector of b1 draws.
                       _["bscale0"]    = bscale0_post,      // nsim-length vector of b0 draws.
                       _["eta_mu"]     = eta_post_con,      // nsim-length vector for eta_mu draws.
                       _["gamma_mu"]   = gamma_post_con,    // nsim-length vector of gamma_mu draws.
                       _["alpha_con"]  = alpha_con,         // matrix of MH alpha's for prognostic trees
                       _["alpha_mod"]  = alpha_mod          // matrix of MH alpha's for trt trees
   ));
}

