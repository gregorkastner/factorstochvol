#' Generic extraction of covariance matrix
#' 
#' Generic function for extracting model-implied covariance matrices, either
#' from the MCMC output, or from the simulated model.
#'
#' @param x An object of class \code{fsvdraws} or \code{fsvsim}.
#' @param ... Arguments to be passed to methods.
#' 
#' @return Structure containing the model-implied covariance matrix.
#'
#' @family generics
#'
#' @seealso covmat.fsvsample covmat.fsvsim
#'
#' @export

covmat <- function(x, ...) {
 UseMethod("covmat")
}


#' Computes the log returns of a vector-valued time series
#'
#' \code{logret} computes the log returns of a multivariate time
#' series, with optional de-meaning.
#'
#' @param dat The raw data, a matrix with \code{n}
#'  (number of timepoints) rows and \code{m}
#'  (number of component series) columns.
#' @param demean Logical value indicating whether the data should
#' be de-meaned.
#' @param standardize Logical value indicating whether the data should
#' be standardized (in the sense that each component series has an empirical
#' variance equal to one).
#' @param ... Ignored.
#'
#' @return Matrix containing the log returns of the (de-meaned)
#' data.
#'
#' @seealso fsvsample logret.matrix logret.data.frame
#'
#' @export

logret <- function(dat, demean = FALSE, standardize = FALSE, ...) {
 UseMethod("logret")
}
