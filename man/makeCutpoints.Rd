\name{makeCutpoints}
\alias{makeCutpoints}
\title{makeCutpoints: Create list of cutpoints for a design matrix X.}
\usage{
makeCutpoints(x)
}
\arguments{
\item{X}{A matrix; the design matrix for the Functional BART model. (Can use the output of \code{\link{makeModelMatrix}}, which expands categorical predictors to columns of indicators when creating design matrix.) Defaults to NULL.}

\item{gridlen}{Length of grid for continuous predictors.  Defaults to 10,000.}
}
\value{
\item{cutpoints}{A list of vectors, where each is a cutpoint for each covariate (column of design matrix X).}
}
\description{
This function creates a list of cutpoints to be used by the funbart fit functions (\code{\link{funbartFit}}, etc.).
}
\examples{
\dontrun{
## Create list of cutpoints.
cuts <- makeCutpoints(X)
}}
\author{
Jennifer E. Starling
}
\keyword{funbart}

