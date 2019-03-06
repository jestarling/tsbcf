#' Simulated cosine data with four covariates.
#'
#' A simulated data set to illustrate the Functional BART
#' method.  Data contains four 'x' covariates; each pair
#' is simulated from bivariate Gaussians with moderate correlation
#' and unit variance.  Each observation has one corresponding time
#' point at which it is observed, randomly chosen from eight discrete
#' time points.
#'
#' Function values are calculated as
#' f(x,t) = g(x1,x2) * cos(t + 2 * pi * h(x2,x3)) where g() and h()
#' are sums of the covariates.
#'
#' This configuration simulates both amplitude change and phase shift for the
#' underlying function.
#'
#' @docType data contained in a csv file with headers
#'
#' @usage read.csv("sim.csv")
#'
#' @format An data frame with n=1000 rows.
#'
#' @keywords datasets
#'
#' @references Starling et al. (2018)
#'
#' @source \href{http://github.com/jestarling/funbart}{funbart R package}
#'
#' @examples
#' data(sim)
#' head(sim)
"sim"
