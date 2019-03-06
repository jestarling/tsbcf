\name{tsbart}
\alias{tsbart}
\title{tsbart: Fits tsbart model for continuous or binary responses.}
\usage{
tsbart(y, tgt, tpred, x, xpred, nburn, nsim, ntree=200,
                   lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                   ecross=1, base_tree=.95, power_tree=2,
                   sd_control=sd(y), use_fscale=TRUE,
                   probit=FALSE, yobs=NULL, verbose=T, mh=F, save_inputs=T)
}
\arguments{
\item{y}{Length n vector with continuous response values.  In probit case, should contain initializations for latent variables.}

\item{tgt}{Length n targeted covariate over which to smooth.}

\item{x}{A (n x p) data frame or matrix containing the covariates which are not to be smoothed over.}

\item{tpred}{Length npred out-of-sample targeted covariate.  Optional.}

\item{xpred}{A (npred x p) data frame or matrix containing out of sample covariates. Optional.}

\item{nburn}{Number of burn-in MCMC iterations. Defaults to 100.}

\item{nsim}{Number of MCMC iterations to save after burn-in. Defaults to 1000.}

\item{ntree}{Number of trees. Defaults to 200.}

\item{lambda}{Scale parameter in the chisq prior sigma^2.  Defaults to NULL, ie being estimated from sigq and sighat. Not appicable for binary case.}

\item{sigq}{Calibration quantile in the chisq prior on sigma^2.  Defaults to .9. Not applicable for binary case.}

\item{sighat}{Calibration estimate for chisq prior on sigma^2.  Defaults to being estimated from data using linear model.  Not applicable for binary case.}

\item{nu}{Degrees of freedom in the chisq prior on sigma^2. Defaults to 3. Not applicable for binary case.}

\item{ecross}{Smoothing parameter; number of expected times f(x,t) crosses the mean response over time, alpha(t). Defaults to 1.  Can set to "tune" to perform parameter tuning.}

\item{base_tree}{Base for tree prior.  Defaults to 0.95.}

\item{power_tree}{Power for tree prior.  Defaults to 2.0.}

\item{sd_control}{SD(f(x,t)) marginally at any covariate value (or its prior median if use_fscale=TRUE).  Default is 2*sd(y).}

\item{use_fscale}{Use a half-Cauchy prior on the scale of f(x,t).}

\item{probit}{0 indicates continuous response; 1 indicates probit response.  Default is 0.  If 1, then yobs must be populated.}

\item{yobs}{Length n vector of binary responses; only populated for probit=1 case.}

\item{verbose}{Boolean for writing progress report to console.}

\item{mh}{Boolean for including Metropolis acceptance detail in output.  Defaults to FALSE.  If true, output includes the metrop dataframe.}

\item{save_inputs}{Boolean for saving user inputs.  If true, output includes the inputs dataframe.}

}
\value{
A list containing the following items:

\item{mcmcdraws}{A (nsim x n) matrix containing the in-sample MCMC draws for the tsbart fit.}

\item{mcmcdraws_oos}{A (nsim x npred) matrix containing the out-of-sample MCMC draws for the tsbart fit, if predicting for out of sample data.}

\item{sigma}{A vector containing the MCMC draws for sigma.}

\item{treefit_sd}{A vector containing draws for the standard deviation of the BART fit f(x,t), equal to sd_control * eta.}

\item{eta}{A vector containing the MCMC draws for eta.}

\item{gamma}{A vector containing the MCMC draws for gamma.}

\item{ecross}{The expected number of crossings.}

\item{metrop}{A ((nburn + nsim)*ntree x 5) dataframe of tree proposal details.  Included if mh=TRUE.  Columns are:
    iter (MCMC iteration number, from 1 to nburn+nsim),
    tree (tree number, from 1 to ntree),
    accept (1=0 accepted, 0=rejected),
    alpha (the Metropolis Hasting alpha for the proposed move), and
    bd ('birth' or 'death', indicating type of tree proposal).)}

\item{inputs}{A dataframe with key function inputs saved.  Excludes data elements.}

}
\description{
This function implements the tsbart method for a continuous or binary responses, smoothing over the covariate tgt,
with unsmoothed covariates x.  Predicts function values at a set of new time points and covariates.
}
\examples{
\dontrun{
See vignette for full example, including expected crossings parameter tuning.
Function call example is as follows.

# Fit model.
fit = tsbart(y, tgt, tpred, x, xpred, nburn, nsim, ntree=200, lambda=NULL, sigq=.9, sighat=NULL,
               nu=3, ecross=1, base_tree=.95, base_tree=.95, power_tree=2,
               sd_control=sd(y), use_fscale=TRUE,
               probit=F, yobs=NULL, verbose=T, mh=F, save_inputs=T)
}}
\author{
Jennifer E. Starling
}
\keyword{tsbart}


