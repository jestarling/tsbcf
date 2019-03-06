\name{survPrep}
\alias{survPrep}
\title{survPrep: Expand a data frame for survival analysis.}
\usage{
survPrep(df, timevar='t', statusvar='status')
}
\arguments{
\item{df}{A data frame to expand for survival analysis.  Must contain a variable with time points > 0, and a binary response column indicating whether the event was observed (1) or not (0) at the observed time. Defaults to NULL.}

\item{timevar}{The name of the column in df containing time values.  Defaults to 't'.}

\item{statusvar}{The name of the column in df containing the binary response. Defaults to 'status'.}
}
\value{
\item{surv}{A data frame containing one row for every time point for every observation, up to the observed time for each observation. The statusvar is set to zero for all times before the observed time.}
}
\description{
This is a utility function for preparing a data frame for survival analysis.
}
\examples{
\dontrun{
   # Test data frame to illustrate.
   testdf = data.frame('id' = c(1,2,3,4,5),
                          't' = c(4,3,1,2,10),
                          'x' = c('a','b','c','d','e'),
                          'y' = c(0,1,1,0,1),
                          'z' = rnorm(5))

   # Create expanded data set.
   survdf = survPrep(testdf, 't', 'y')
}}
\author{
Jennifer E. Starling
}
\keyword{funbart, dataprep}

