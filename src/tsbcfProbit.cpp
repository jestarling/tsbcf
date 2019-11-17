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
List tsbcfProbit(arma::vec y,       // True latent variable values of the response.
              arma::vec yobs,       // Vector of observed y values; 1 and 0.
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
              int nburn,
              int nsim,
              int ntree_con=200,
              int ntree_mod=50,
              double lambda=-999,
              double sigq=.9,
              double nu=3,
              double offset=0,
              double base_con=.95,
              double power_con=2.0,
              double base_mod=.25,
              double power_mod=3.0,
              double ecross_con=1,
              double ecross_mod=1,
              double con_sd=1,
              double mod_sd=1,
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
   // Read, format info about tgt
   //*****************************************************************************

   vec tref = unique(join_cols(tgt,tpred));  // Vector of unique tgt points, sorted in asc. order.
   int tlen = tref.n_elem;                   // Number of unique tgt points (T).

   if(silent_mode==false){
      Rcout << "unique tgt points: " << endl << tref << endl;
      Rcout << "Number of total tgt points T: " << tlen << endl;
   }

   //------------------------------------------------------
   // Calculate alpha_t (ybar_t) vector, for conditional-on-tgt mean-centering data.
   //------------------------------------------------------

   vec alpha_t = zeros(tlen);    // Holds ybar at each t.
   vec n_t = zeros(tlen);        // Holds sample sizes at each t.
   uvec idx;                     // Holds indices where ti = t, for calculating ybar_t.

   // Iterate through tgt points.
   for(int i = 0; i < tlen; i++){

      // All observations.
      idx = find(tgt == tref(i));
      n_t(i) = idx.size();

      if(n_t(i) > 0){
         alpha_t(i) = sum(y(idx)) / n_t(i);
      }
   }

   // Create n-length vector with corresponding ybar_i for each y_i.
   // Will let us subtract off tgt-means and add back later.
   vec alpha_ti(y.size());
   idx.reset();

   for(size_t i = 0; i < y.size(); i++){
      idx = find(tref == tgt(i));
      alpha_ti(i) = as_scalar(alpha_t(idx));
   }

   //*****************************************************************************
   // Format y, and calculate sufficient statistics for y, for initializing
   //*****************************************************************************

   // Mean-center data, conditional on the mean for each tgt point.
   y = y - alpha_ti;

   // Scalar versions of sufficient stats for all y.
   sinfo allys;
   allys.n = n;
   allys.sy = sum(y);
   allys.sy2 = sum(y % y);

   // Vector sufficient stats for y_t's (T-length tgt vector).
   vec ybar_vec = zeros(tlen);
   allys.n_vec = n_t;
   allys.sy_vec = zeros(tlen);

   // Iterate through tgt points; calculate ybar_t.
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

   // Number of in-sample control, treatment covariates.
   size_t p_con = x_con.size()/n;
   size_t p_mod = x_mod.size()/n;

   // Number of out-of-sample control, treatment covariates.
   size_t p_con_pred = xpred_con.size()/n_pred;
   size_t p_mod_pred = xpred_mod.size()/n_pred;

   // Check correct number of control, treatment covariates.
   double npred_con = xpred_con.size()/p_con;
   double npred_mod = xpred_mod.size()/p_mod;
   if(xpred_con.size() != npred_con*p_con) Rcout << "error, wrong number of elements in prediction data set\n";
   if(xpred_mod.size() != npred_mod*p_mod) Rcout << "error, wrong number of elements in prediction data set for treatment effects\n";

   //------------------------------------------------------------------------
   //x cutpoints
   //------------------------------------------------------------------------

   xinfo xi_con;
   xinfo xi_mod;

   // control
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
   // Initialize tree fit scale-related variables.  These will be updated as we go.
   //------------------------------------------------------------------------

   double bscale_prec = 2;
   double bscale0 = -0.5;
   double bscale1 = 0.5;

   double mscale_prec = 1.0;
   double mscale = 1.0;
   double delta_con = 1.0;
   double delta_mod = 1.0;

   //------------------------------------------------------------------------
   // Initialize trees.
   //------------------------------------------------------------------------

   // mu(x,t,pi) for control trees. if you sum the fit over the trees you get the fit.
   std::vector<tree> t_con(ntree_con);
   for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar_vec / ntree_con);

   // tau(x,t) for treatment trees
   std::vector<tree> t_mod(ntree_mod);
   for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init / ntree_mod);

   //------------------------------------------------------------------------
   // Prior for control trees.
   //------------------------------------------------------------------------
   pinfo pi_con(tlen);

   pi_con.pbd = 1.0;          //prob of birth/death move
   pi_con.pb = .5;            //prob of birth given  birth/death
   pi_con.alpha = base_con;   //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi_con.beta = power_con;   //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node

   // Sigma.  Use user-provided sigma, unless null (-999); then use sd(y).  Gets adjusted later for backfitting: sig^2 / eta_con^2.
   double sigma = 1;
   pi_con.sigma = 1;
   pi_con.sigma = sigma / fabs(mscale); //resid variance in backfitting is \sigma^2_y/mscale^2

   // Use override for lambda if provided by user.  Else calculate using sigma and sigq.
   double qchi = 0;
   if(lambda==-999){
      qchi = R::qchisq(1-sigq, nu, true, false);
      lambda = (sigma * sigma * qchi) / nu;
   }

   // Initialize new tgt-related prior components.
   pi_con.mu0 = zeros(tlen);
   pi_con.ls = tlen / (PI * ecross_con);

   // Initialize GP-related prior components.
   pi_con.mu0 = zeros(tlen);
   pi_con.ls = tlen / (PI * ecross_con);

   pi_con.tau = con_sd / sqrt(delta_con) * sqrt((double) ntree_con);
   pi_con.K = cov_se(tref, tref, pi_con.ls, 1).i();
   pi_con.Prec = (1/(pi_con.tau*pi_con.tau)) * pi_con.K; // See info.h line 67.

   //------------------------------------------------------------------------
   // Prior for treatment trees.
   //------------------------------------------------------------------------
   pinfo pi_mod(tlen);

   pi_mod.pbd = 1.0;                // prob of birth/death move
   pi_mod.pb = .5;                  // prob of birth given  birth/death
   pi_mod.alpha = base_mod;         // prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi_mod.beta = power_mod;         // 3 for treatment trees.
   pi_mod.sigma = pi_con.sigma;     // Same initial value as pi_con.sigma.  Gets adjusted later for backfitting: sig^2 / eta_mod^2.

   // Initialize new tgt-related prior components.
   pi_mod.mu0 = zeros(tlen);
   pi_mod.ls = tlen / (PI * ecross_mod);

   // Initialize GP-related prior components.
   pi_mod.mu0 = zeros(tlen);
   pi_mod.ls = tlen / (PI * ecross_mod);

   pi_mod.tau = con_sd / sqrt(delta_mod) * sqrt((double) ntree_mod);
   pi_mod.K = cov_se(tref, tref, pi_mod.ls, 1).i();
   pi_mod.Prec = (1/(pi_mod.tau*pi_mod.tau)) * pi_mod.K; // See info.h line 67.


   //------------------------------------------------------------------------
   // dinfo for control function mu(x,t)
   //------------------------------------------------------------------------
   double* allfit_con = new double[n];    // Holds sum of fit of all trees for each observation.
   double* r_con = new double[n];         // Holds residuals.  y - (allfit - ftemp) = y - allfit + ftemp.

   // Initialize allfit_con to the ybar_t for each obs' corresponding tgt value.
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

   // Initialize allfit to the trt_init_t = (ybar_trt_t - ybar_ctrl_t) for each obs' corresponding tgt value.
   for (int i=0; i<n; ++i){
      idx = find(tref==tgt(i));
      allfit_mod[i] = as_scalar( (z[i]*bscale1 + (1-z[i])*bscale0) * trt_init(idx));
   }

   dinfo di_mod;
   di_mod.n=n; di_mod.p=p_mod; di_mod.x = &x_mod[0]; di_mod.y=r_mod; //the "y" for each draw will be the residual
   di_mod.tlen = tlen; di_mod.tref = tref; di_mod.t = &tgt[0];

   //------------------------------------------------------------------------
   // dinfo for control function mu(x,t) out-of-sample
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

   // Standard deviations for mu and tau.
   NumericVector mu_sd_post(nsim);
   NumericVector tau_sd_post(nsim);

   // For control scale.
   NumericVector mscale_post(nsim);

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

   // For holding MH alphas for control and treatment trees.
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
      //draw control trees
      //------------------------------------------------------------------------

      // Loop through ntree_con control trees.
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
            allfit[k] = allfit[k] - mscale*ftemp[k];
            allfit_con[k] = allfit_con[k] - mscale*ftemp[k];
            r_con[k] = (y[k]-allfit[k])/mscale;

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
            allfit[k] += mscale*ftemp[k];
            allfit_con[k] += mscale*ftemp[k];
         }

      } // End control tree loop.


      //------------------------------------------------------------------------
      //draw treatment trees
      //------------------------------------------------------------------------

      // Set up temporary 'precisions' for heteroskedastic updates.
      for(size_t k=0;k<ntrt;k++){
         precs[k] = (bscale1*bscale1)/(sigma*sigma);
      }
      for(size_t k=ntrt;k<n;k++){
         precs[k] = (bscale0*bscale0)/(sigma*sigma);
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
      // For control: Update mscale and delta_con.
      //------------------------------------------------------------------------

      if(use_muscale==true){

         double ww = 0.0;  // Tracks sum of w*w
         double rw = 0.0;  // Tracks sum of r*w
         double s2 = sigma*sigma;

         for(size_t k=0; k<n; ++k) {
            double w = s2*mscale*mscale/(allfit_con[k]*allfit_con[k]);
            if(w!=w) {
               Rcout << " w " << w << endl;
               stop("");
            }

            double r = (y[k] - allfit_mod[k])*mscale/allfit_con[k];
            if(r!=r) {
               stop("");
            }
            ww += 1/w;
            rw += r/w;
         }

         double mscale_old = mscale;
         double mscale_fc_var = 1/(ww + mscale_prec);
         mscale = mscale_fc_var*rw + gen.normal(0., 1.)*sqrt(mscale_fc_var);

         for(size_t k=0; k<n; ++k) {
            allfit_con[k] = allfit_con[k]*mscale/mscale_old;
         }

         // update delta_con
         double ssq = 0.0;
         tree::npv bnv;
         typedef tree::npv::size_type bvsz;
         double endnode_count = 0.0;

         for(size_t j=0;j<ntree_con;j++) {
            bnv.clear();
            t_con[j].getbots(bnv);
            bvsz nb = bnv.size();
            for(bvsz ii = 0; ii<nb; ++ii) {
               arma::vec mm = bnv[ii]->getm(); //node parameter
               ssq += as_scalar(mm.t() * pi_con.Prec * mm);
               //ssq += mm*mm/(pi_con.tau*pi_con.tau);
               endnode_count += 1.0;
            }
         }

         //delta_con = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
         delta_con = gen.gamma(0.5*(15.0 + tlen*endnode_count), 1.0)/(0.5*(15.0 + ssq));
         pi_con.tau = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con));
         pi_con.Prec = (1 / (pi_con.tau * pi_con.tau)) * pi_con.K;

      } else {
         mscale = 1.0;
      }

      //------------------------------------------------------------------------
      // Update bscale.
      //------------------------------------------------------------------------

      if(use_tauscale==true){

         // Calculate sums required for means and variances.
         double ww0 = 0.0, ww1 = 0.0; // holds precisions.
         double rw0 = 0.0, rw1 = 0.0; // holds resids*precisions.
         double s2 = sigma * sigma;

         for(size_t k=0;k<n;k++){
            double bscale = (k<ntrt) ? bscale1 : bscale0;
            double w = s2*bscale*bscale / (allfit_mod[k]*allfit_mod[k]);

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
               ww1 += 1/w;
               rw1 += r/w;
            } else{
               ww0 += 1/w;
               rw0 += r/w;
            }
         }

         // Calculate mean and variance for b0 and b1 full conditionals, and draw.
         double bscale1_old = bscale1;
         double bscale_fc_var = 1/(ww1 + bscale_prec);
         bscale1 = bscale_fc_var*rw1 + gen.normal(0., 1.)*sqrt(bscale_fc_var);

         double bscale0_old = bscale0;
         bscale_fc_var = 1/(ww0 + bscale_prec);
         bscale0 = bscale_fc_var*rw0 + gen.normal(0., 1.)*sqrt(bscale_fc_var);

         // Rescale fits.
         for(size_t k=0;k<ntrt;k++){
            allfit_mod[k] = allfit_mod[k]*bscale1/bscale1_old;
         }
         for(size_t k=ntrt;k<n;k++){
            allfit_mod[k] = allfit_mod[k]*bscale0/bscale0_old;
         }

         // update delta_mod
         double ssq = 0.0;
         tree::npv bnv;
         typedef tree::npv::size_type bvsz;
         double endnode_count = 0.0;
         double df_b = 30;  // Df for delta ~ Ga(df/2, df/2) prior.

         for(size_t j=0;j<ntree_mod;j++) {
            bnv.clear();
            t_mod[j].getbots(bnv);
            bvsz nb = bnv.size();
            for(bvsz ii = 0; ii<nb; ++ii) {
               arma::vec mm = bnv[ii]->getm(); //node parameter
               ssq += as_scalar(mm.t() * pi_con.Prec * mm);
               endnode_count += 1.0;
            }
         }

         //delta_mod = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
         delta_mod = gen.gamma(0.5*(df_b + tlen*endnode_count), 1.0)/(0.5*(df_b + ssq));
         pi_mod.tau   = mod_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod));
         pi_mod.Prec = (1 / (pi_mod.tau * pi_mod.tau)) * pi_mod.K;

      } else{
         bscale1 = .5;
         bscale0 = -.5;
      }

      pi_mod.sigma = sigma;

      //------------------------------------------------------------------------
      // Sync yhat=allfit[k] after scale updates.
      //------------------------------------------------------------------------
      if(use_muscale || use_tauscale){
         for(size_t k=0;k<n;k++){
            allfit[k] = allfit_con[k] + allfit_mod[k];
         }
      }

      //------------------------------------------------------------------------
      // Draw sigma.  No draws here since probit, but still have to update scale for backfitting.
      //------------------------------------------------------------------------

      pi_con.sigma = sigma / fabs(mscale);
      pi_mod.sigma = sigma;

      //------------------------------------------------------------------------
      // Probit draw.
      //------------------------------------------------------------------------
      for(size_t k=0; k<n; ++k) {
         if(yobs[k]>0) {
            y[k] = rtnormlo1(allfit[k], -offset); //offset isn't in allfit
         } else {
            y[k] = rtnormhi1(allfit[k], -offset);
         }
      }

      //------------------------------------------------------------------------
      // Save MCMC output.
      //------------------------------------------------------------------------

      // First, need a vector for the alpha_ti values at each tpred.
      // For out of sample observations, pull in appropriate ybar_t
      // from each tgt point.  (Uses ybar_t from in-sample data for each tgt point).
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

         // Save traces of important parameters.
         sigma_post(i-nburn) = 1.0;
         mscale_post(i-nburn) = mscale;
         bscale1_post(i-nburn) = bscale1;
         bscale0_post(i-nburn) = bscale0;
         mu_sd_post(i-nburn) = fabs(mscale)*con_sd;
         tau_sd_post(i-nburn) = fabs(bscale1-bscale0)*mod_sd;

         // Save in-sample treatment effects, control effects, and yhat fits.
         for(size_t k=0;k<n;k++) {
            yhat_post(i-nburn,k) = allfit[k] + alpha_ti[k];
            double bscale = (k<ntrt) ? bscale1 : bscale0;
            mod_post(i-nburn,k) = (bscale1-bscale0) * allfit_mod[k]/bscale; // allfit_mod[k] = [b1*z_i + b0*(1-z_i)] * tau(x,t) == bscale * tau(x,t) --> need (b1 - b0) * tau(x,t), so divide out bscale.
            con_post(i-nburn,k) = yhat_post(i-nburn,k) - z[k]*mod_post(i-nburn,k);
         }

         // Draw out of sample treatment effects.
         if(di_mod_pred.n>0){
            for(size_t k=0;k<n_pred;k++) {
               double bscale = (k<ntrt_pred) ? bscale1 : bscale0;
               mod_est_post(i-nburn,k) = (bscale1 - bscale0) *fit_i(k,t_mod,xi_mod,di_mod_pred);

               yhat_est_post(i-nburn,k) = alpha_ti_pred[k] +
                  mscale*fit_i(k,t_con,xi_con,di_con_pred) +
                  bscale*fit_i(k,t_mod,xi_mod,di_mod_pred);
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
                       _["mu"]         = con_post,          // (nsim x n) matrix of in-sample mcmc control effects.
                       _["mu_oos"]     = con_est_post,      // (nsim x npred) matrix of out-of-sample mcmc control effects.
                       _["tau"]        = mod_post,          // (nsim x n) matrix of in-sample mcmc treatment effects.
                       _["tau_oos"]    = mod_est_post,      // (nsim x npred) matrix of out-of-sample mcmc treatment effects.
                       _["sigma"]      = sigma_post,        // nsim-length vector of sigma draws.
                       _["mscale"]     = mscale_post,       // nsim-length vector of mscale draws.
                       _["bscale1"]    = bscale1_post,      // nsim-length vector of b1 draws.
                       _["bscale0"]    = bscale0_post,      // nsim-length vector of b0 draws.
                       _["mu_sd_post"] = mu_sd_post,        // Marginal sd of control tree fit.
                       _["tau_sd_post"]= tau_sd_post,       // Marginal sd of treatment tree fit.
                       _["alpha_con"]  = alpha_con,         // matrix of MH alpha's for control trees
                       _["alpha_mod"]  = alpha_mod ,         // matrix of MH alpha's for trt trees
                       _["alpha_t"]   = alpha_t
   ));
}

