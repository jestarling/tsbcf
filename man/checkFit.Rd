\name{checkFit}
\alias{checkFit}
\title{checkFit: Generates log-likelihood and WAIC for tsbart fits.  Called by tuneEcross function.}
\usage{
checkFit(y, mcmcdraws, probit=FALSE, doWaic=TRUE, sig, yobs=NULL)
}
\arguments{
\item{y}{A vector with continuous response values. Required, defaults to NULL.}

\item{mcmcdraws}{A matrix containing the MCMC draws for the tsbart fit. Rows are MCMC iterations, columns are observations.  Number of columns must match length of y.  Required, defaults to NULL.}


\item{probit}{Boolean indicator for whether model is probit. Optional, defaults to FALSE.}

\item{doWaic}{Boolean indicator for whether to calculate WAIC in addition to log-likelihood. Optional, defaults to TRUE.}

\item{sig}{A vector with sigma MCMC draws.  Required if probit=FALSE, defaults to NULL.  Must be NULL if probit=TRUE.}

\item{yobs}{A vector of 0/1 indicators for each observation.  Required if probit=TRUE, defaults to NULL.}
}
\value{
A list containing the following components:

\item{waic}{The waic value for the model fit.}

\item{ll_mcmc}{The log-likelihood for all mcmc draws.}
}
\description{
This function calculates the log-likelihood of a tsbart model fit, and optionally the WAIC.
}
\examples{
\dontrun{
# Function call example is as follows.  Let fit be the output of tsbart().

## In-sample example:
  checkFit(y=y,
         mcmcdraws = fit[["mcmcdraws"]],
         probit=FALSE, doWaic=TRUE,
         sig = fit[["sigma"]],
         yobs = NULL)

## Out-of-sample example:
  checkFit(y=y_pred,
         mcmcdraws = fit[["mcmcdraws_oos"]],
         probit=FALSE, doWaic=TRUE,
         sig = fit[["sigma"]],
         yobs = NULL)
}

}
\author{
Jennifer E. Starling
}
\keyword{tsbart}


