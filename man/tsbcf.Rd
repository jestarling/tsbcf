\name{tsbcf}
\alias{tsbcf}
\title{tsbcf: Fits the Targeted Smooth Bayesian Causal Forests model for continuous or binary
responses.}
\usage{
tsbcf(y, pihat, z, tgt, x_control, x_moderate,
                  pihatpred=NULL, zpred=NULL, tpred=NULL, xpred_control=NULL, xpred_moderate=NULL,
                  nburn=100, nsim=1000, ntree_control=200, ntree_moderate=50,
                  lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                  base_control=.95, power_control=2,
                  base_moderate=.25, power_moderate=3,
                  sd_control=2*sd(y), sd_moderate=sd(y),
                  treatment_init = rep(1,length(unique(tgt))),
                  use_muscale=T, use_tauscale=T,
                  ecross_control=1, ecross_moderate=1,
                  ecross_control_candidates = NULL, ecross_moderate_candidates = NULL,
                  ecross_tune_nsim=100, ecross_tune_nburn=1000,
                  pihat_in_trt=F,
                  link='identity', yobs=NULL, set_scales=F,
                  verbose=T, mh=F, save_inputs=T)
}
\arguments{
\item{y}{Length n vector with continuous response values.  In probit case, should contain
initializations for latent variables (+/- 1.96, for example).}

\item{pihat}{Length n vector with estimated propensity score values.}

\item{z}{Length n vector with binary treatment indicators (must be 1/0).}

\item{tgt}{Length n targeted covariate over which to smooth.}

\item{x_control}{A (n x p_control) data frame or matrix containing the covariates which are not to
be smoothed over, for modeling prognostic effect.}

\item{x_moderate}{A (n x p_moderate) data frame or matrix containing the covariates which are not
to be smoothed over, for modeling treatment effect.}

\item{pihatpred}{Length npred vector with out-of-sample estimated propensity score values.}

\item{zpred}{Length npred vector with out-of-sample binary treatment indicators (must be 1/0).}

\item{tpred}{Length npred out-of-sample targeted covariate.}

\item{xpred_control}{A (npred x p_control) data frame or matrix containing out of sample covariates,
for modeling prognostic effect.}

\item{xpred_moderate}{A (npred x p_moderate) data frame or matrix containing out of sample covariates,
for modeling treatment effect.}

\item{nburn}{Number of burn-in MCMC iterations. Defaults to 100.}

\item{nsim}{Number of MCMC iterations to save after burn-in. Defaults to 1000.}

\item{ntree_control}{Number of trees for prognostic fit. Defaults to 200.}

\item{ntree_moderate}{Number of trees for treatment fit. Defaults to 50.}

\item{lambda}{Scale parameter in the chisq prior sigma^2.  Defaults to NULL, ie being estimated
from sigq and sighat. Not appicable for binary case.}

\item{sigq}{Calibration quantile in the chisq prior on sigma^2.  Defaults to .9. Not applicable for
binary case.}

\item{sighat}{Calibration estimate for chisq prior on sigma^2.  Defaults to being estimated from
data using linear model.  Not applicable for binary case.}

\item{nu}{Degrees of freedom in the chisq prior on sigma^2. Defaults to 3. Not applicable for
binary case.}

\item{base_control}{Base for prognostic tree prior.  Defaults to 0.95.}

\item{power_control}{Power for prognostic tree prior.  Defaults to 2.}

\item{base_moderate}{Base for treatment tree prior.  Defaults to 0.25.}

\item{power_moderate}{Power for treatment tree prior.  Defaults to 3.}

\item{sd_control}{Marginal scale for the mu(x,t) control fit at any covariate value.  Default is
2*sd(y).  If use_muscale=F, corresponds to sigma_mu in original BART model.}

\item{sd_moderate}{Marginal scale for the tau(x,t) treatment fit at any covariate value.  Default is
2*sd(y). If use_tauscale=F, corresponds to sigma_mu in original BART model.}

\item{treatment_init}{T-length vector of initial treatment values.  Defaults to 1's.  Must be
length of number of unique values in tgt and tpred combined.}

\item{use_muscale}{Use a half-Cauchy prior on the scale of mu(x,t), the prognostic fit.}

\item{use_tauscale}{Use a half-Normal prior on the scale of tau(x,t), the treatment fit.}

\item{ecross_control}{Smoothing parameter for prognostic tree fit; number of expected times f(x,t)
crosses the mean response over time, alpha(t). Defaults to 1.  Can set to "tune" to perform parameter tuning.}

\item{ecross_moderate}{Smoothing parameter for treatment tree fit; number of expected times f(x,t)
crosses the mean response over time, alpha(t). Defaults to 1.  Can set to "tune" to perform parameter tuning.}

\item{ecross_control_candidates}{Vector of candidate values for ecross_control.  Defaults to NULL;
is only considered if ecross_control="tune", otherwise is ignored.}

\item{ecross_moderate_candidates}{Vector of candidate values for ecross_moderate.  Defaults to NULL;
is only considered if ecross_moderate="tune", otherwise is ignored.}

\item{ecross_tune_nburn}{Number of burn-in MCMC iterations for tuning the expected crossing
parameters. Defaults to 100.}

\item{ecross_tune_nsim}{Number of MCMC iterations for tuning the expected crossing parameters.
Defaults to 1000.}

\item{pihat_in_trt}{Boolean, ndicates whether to include the estimated propensity score as a
covariate for the treatment tree fit.  (Is always included in prognostic tree fit.)  Defaults to TRUE.}

\item{link}{Set to "identity" for Gaussian response, "probit", or "logit".  Default is "identity".  If not "identity",
then yobs must be populated.}

\item{yobs}{Length n vector of binary responses; only populated for probit=1 case.}

\item{set_scales}{Defaults to F.  When probit=T, if true, then sets control_sd and
moderate_sd based on estimates of baseline and relative risk from data.}

\item{verbose}{Boolean for writing progress report to console.}

\item{mh}{Boolean for including Metropolis acceptance detail in output.  Defaults to FALSE.  If
true, output includes the metrop_control and metrop_moderate dataframes.}

\item{save_inputs}{Boolean for saving user inputs.  If true, output includes the inputs dataframe.}

}
\value{
A list containing the following items:

\item{tgt}{A n-length vector containing the user input for tgt.}

\item{yhat}{A (nsim x n) matrix containing the in-sample MCMC draws for the tsbcf fit.}

\item{mu}{A (nsim x n) matrix containing the in-sample MCMC draws for the prognostic portion of the
tsbcf fit.}

\item{tau}{A (nsim x n) matrix containing the in-sample MCMC draws for the treatment portion of the
tsbcf fit.}

\item{tpred}{A pred-length vector containing the user input for tpred.}

\item{yhat_oos}{A (nsim x npred) matrix containing the out-of-sample MCMC draws for the tsbcf fit.}

\item{mu_oos}{A (nsim x npred) matrix containing the out-of-sample MCMC draws for the prognostic
portion of the tsbcf fit.}

\item{tau_oos}{A (nsim x npred) matrix containing the out-of-sample MCMC draws for the treatment
portion of the tsbcf fit.}

\item{sigma}{A nsim-length vector containing the MCMC draws for sigma.}

\item{mu_sd}{A nsim-length vector containing the MCMC draws for the standard deviation of the
prognostic tree fits.}

\item{tau_sd}{A nsim-length vector containing the MCMC draws for the standard deviation of the
treatment tree fits.}

\item{ecross_control}{The expected number of crossings for the prognostic tree fit.}

\item{ecross_moderate}{The expected number of crossings for the treatment tree fit.}

\item{metrop}{A ([(nburn + nsim)*ntree_control + (nburn + nsim)*ntree_moderatel] x 6) dataframe of
tree proposal details.  Included if mh=TRUE.  Columns are:
    iter (MCMC iteration number, from 1 to nburn+nsim),
    tree (tree number, from 1 to ntree),
    accept (1=0 accepted, 0=rejected),
    alpha (the Metropolis Hasting alpha for the proposed move),
    bd ('birth' or 'death', indicating type of tree proposal), and
    tree ('control' or 'moderate').}

\item{inputs}{A dataframe with key function inputs saved.  Excludes data elements.}

}
\description{
This function implements the tsbcf method for a continuous or binary responses, smoothing over
the covariate tgt,
with unsmoothed covariates x_control for prognostic tree fit and x_moderate for treatment tree
fit.  Predicts function values at a set of new time points and covariates.
}
\examples{
\dontrun{
# Function call example is as follows.
fit = tsbcf(y, pihat, z, tgt, x_control, x_moderate,
                  pihatpred=NULL, zpred=NULL, tpred=NULL, xpred_control=NULL, xpred_moderate=NULL,
                  nburn=100, nsim=1000, ntree_control=200, ntree_moderate=50,
                  lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                  base_control=.95, power_control=2,
                  base_moderate=.25, power_moderate=3,
                  sd_control=2*sd(y), sd_moderate=sd(y),
                  treatment_init = rep(1,length(unique(tgt))),
                  use_muscale=T, use_tauscale=T,
                  ecross=1, ecross_mod=1,
                  ecross_control_candidates = NULL, ecross_moderate_candidates = NULL,
                  ecross_tune_nsim=100, ecross_tune_nburn=1000,
                  pihat_in_trt=T,
                  link='identity', yobs=NULL, set_scales=F,
                  verbose=T, mh=F, save_inputs=T)
}}
\author{
Jennifer E. Starling
}
\keyword{tsbcf}


