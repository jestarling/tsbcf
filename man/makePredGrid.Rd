\name{makePredGrid}
\alias{makePredGrid}
\title{makePredGrid: Expand a data frame of predicted observations over a grid.}
\usage{
makePredGrid(df_new, timevar, timegrid)

  #           timevar  = name of the time variable, for formatting output data set.
   #           timegrid = grid of timepoints at which to predict.
}
\arguments{
\item{df_new}{A data frame of observations at which to predict. Defaults to NULL.}

\item{timevar}{The name of the variable over which to grid (such as time).  In quotes, e.g. timevar = "t".  Defaults to NULL.}

\item{timegrid}{A vetor of timepoints over which to grid.  Defaults to NULL.}
}
\value{
\item{predgrid}{A data frame containing one row for each observation in df_new, at each unique grid value of timevar included in timegrid.  Includes an "obs_id" column to identify each unique observation, which has one row per grid value.}
}
\description{
This is a utility function for taking a data frame of new observations at which to predict, and expanding the data frame to include one observation at each unique grid value.  Example: Input a data frame of patient characteristics, and get a data frame with one row for each patient at each point in time.

For predicting at new observations in funbartFit and funbartProbit, the user is not required to use output from this function; it simply provides a quick way to expand observations over some grid.
}
\examples{
\dontrun{
## Create grid
new_obs_grid <- makePredGrid(df_new, "t", timegrid=seq(1,10,by=1))
}}
\author{
Jennifer E. Starling
}
\keyword{funbart}

