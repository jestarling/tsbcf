\name{getCausalEstimands}
\alias{getCausalEstimands}
\title{getCausalEstimands}
\usage{
getCausalEstimands(tsbcf_output, subgroups=NULL, probit=F, relrisk=F)
}
\arguments{
\item{tsbcf_obj}{Output from the tsbcf() function, run with probit=T.}

\item{subgroups}{An n-length vector of subgroup values; these can be numeric values, factors, or strings.  Defaults to NULL.  If populated, CATE and CATE_t are returned instead of ATE and ATE_t.}


\item{probit}{Indicator for whether the tsbcf() output was generated as a probit model.  Defaults to FALSE.}

\item{relrisk}{Indicator, for probit models, of whether to return the ATE or the relative risk.  Note, for probit models, these are calculated on probability scale; w1 = pnorm(mu + tau), w0 = pnorm(mu), ATE = w1-w0, RR = w1/w0.  Defaults to FALSE; is ignored if probit=F.}

\item{indiv}{Indicator for whether to include individual treatment effect estimates.  Defaults to FALSE.}
}
\value{

If subgroups = NULL, a list containing the following components:

\item{ate_post}{A (nsim x n) matrix of posterior average treatment effect draws, where (nsim x n) correspond to the tsbcf_output object.}

\item{ate_t_post}{A (nsim x tlen) matrix of posterior average treatment effect draws, where (nsim x tlen) correspond to the tsbcf_output object.  Each row is an MCMC iteration; each column is a unique tgt point.  This is a matrix of posterior draws averaged over observations within each tgt value.}

\item{ate_t_hat}{A (tlen x 4) dataframe of posterior summary information for ate_t.  Columns are: tgt, ate=posterior mean, lb=posterior credible interval lower bound, ub=posterior credible interval upper bound.  Rows are tgt values.  Example: ate=apply(ate_t_post,2,mean), and lb=apply(ate_t_post,2,function(x) quantile(x,.025))}

If subgroups!=NULL, a list containing the following components:

\item{ate_post}{A (nsim x ngroups) matrix of posterior conditional average treatment effect draws, where (nsim x ngroups) correspond to the tsbcf_output object and the number of unique subgroups provided.  Each row is an MCMC iteration; each column is a subgroup.}

\item{ate_t_post}{A (nsim x tlen*ngroups) matrix of posterior conditional average treatment effect draws, where (nsim x tlen*ngroups) correspond to the tsbcf_output object and the number of unique subgroups provided.  Each row is an MCMC iteration; each column is a subgroup-tgt value combination.}

\item{ate_t_hat}{A (tlen*ngroups x 5) dataframe of posterior summary information for cate_t.  Columns are: tgt, subgroup, ate=posterior mean, lb=posterior credible interval lower bound, ub=posterior credible interval upper bound.  Rows are tgt-subgroup combinations. Example: ate=apply(ate_t_post,2,mean), and lb=apply(ate_t_post,2,function(x) quantile(x,.025))}

If indiv=TRUE, the following items are also included:

\item{ate_indiv}{A data frame with n rows and 4 columns, containing individual posterior treatment effect estimates.  Columns are: obs (observation number), posterior mean, lower bound (lb), and upper bound (ub) for the posterior credible interval.  Included if indiv=TRUE.}

\item{ate_indiv_oos}{A data frame with npred rows and 4 columns, containing out-of-sample individual posterior treatment effect estimates.  Columns are: obs (observation number), posterior mean, lower bound (lb), and upper bound (ub) for the posterior credible interval.  Included if indiv=TRUE and out-of-sample data is provided for prediction.}

}
\description{
This function returns the posterior draws, and posterior summary information at each tgt value (posterior mean and credible interval bounds), for the ATE.  If subgroups are given, posterior summary information is given at each tgt value/subgroup combination, for the CATE.  Optionally, individual posterior information is also provided.
}
\examples{
\dontrun{
# Function call example is as follows.  Let fit be the output of tsbcf(), run with probit=T.
  getCausalEstimands(fit, subgroups=NULL, probit=T, relrisk=F)
}

}
\author{
Jennifer E. Starling
}
\keyword{tsbart}


