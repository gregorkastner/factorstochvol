#' Generic extraction of covariance matrix
#' 
#' Generic function for extracting model-implied covariance matrices, either
#' from the MCMC output, or from the simulated model. Details about the
#' function's behavior can be found in  \code{\link{covmat.fsvdraws}}
#' (the function invoked when applied to MCMC output) or
#' \code{\link{covmat.fsvsim}} (the function invoked when applied to a
#' simulated model.
#'
#' @param x An object of class \code{fsvdraws} or \code{fsvsim}.
#' @param ... Arguments to be passed to methods.
#' 
#' @return Structure containing the model-implied covariance matrix.
#'
#' @family generics
#'
#' @export

covmat <- function(x, ...) {
 UseMethod("covmat")
}


#' Generic extraction of correlation matrix
#' 
#' Generic function for extracting model-implied correlation matrices, either
#' from the MCMC output, or from the simulated model. Details about the
#' function's behavior can be found in  \code{\link{cormat.fsvdraws}}
#' (the function invoked when applied to MCMC output) or
#' \code{\link{cormat.fsvsim}} (the function invoked when applied to a
#' simulated model.
#'
#' @param x An object of class \code{fsvdraws} or \code{fsvsim}.
#' @param ... Arguments to be passed to methods.
#' 
#' @return Structure containing the model-implied covariance matrix.
#'
#' @family generics
#'
#' @export

cormat <- function(x, ...) {
 UseMethod("cormat")
}
