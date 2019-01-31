#' Extract "true" model-implied covariance matrix for several points in time
#'
#' \code{covmat} extracts the model-implied (time-varying) covariance matrix
#' from an \code{fsvsim} object.
#'
#' @param x Object of class \code{'fsvsim'}, usually resulting from a call
#' of the function \code{\link{fsvsim}}.
#' @param these Vector indicating which points in time should be extracted,
#' defaults to all.
#' @param ... Ignored.
#'
#' @note Currently crudely implemented as an R loop over all time points,
#' may be slow.
#'
#' @return Array of dimension \code{m} times \code{m} times
#' \code{length(these)}, containing the model-implied covariance matrix.
#'
#' @family simulation
#'
#' @seealso covmat covmat.fsvsample
#'
#' @export

covmat.fsvsim <- function(x, these = seq_len(nrow(x$y)), ...) {
 if (!is(x, "fsvsim")) stop("Argument 'x' must be of class 'fsvsim'.")
 
 if (!is.numeric(these) || min(these) < 1 || max(these) > nrow(x$y))
  stop("Illegal argument value 'these'.")
 
 m <- nrow(x$facload)
 r <- ncol(x$facload)
 covmat <- array(NA_real_, dim = c(m, m, length(these)))
 facload <- x$facload
 
 for (j in seq_along(these)) {
  facvar <- exp(x$facvol[these[j],])
  idivar <- exp(x$idivol[these[j],])
  covmat[,,j] <- tcrossprod(sweep(facload, 2, facvar, '*'), facload)
  diag(covmat[,,j]) <- diag(covmat[,,j]) + idivar
 }

 covmat
}


#' Extract "true" model-implied covariances of two series only
#'
#' \code{covelement} extracts the model-implied (time-varying) covariances between
#' (exactly) two component series.
#'
#' @param x Object of class \code{'fsvsim'}, usually resulting from a call
#' of the function \code{\link{fsvsim}}.
#' @param i Index of component series 1.
#' @param j Index of component series 2.
#' @param these Vector indicating which points in time should be extracted,
#' defaults to all.
#'
#' @return Vector with the requested covariances.
#'
#' @family simulation
#'
#' @export

covelement <- function(x, i, j, these = seq_len(nrow(x$y))) {
 if (!is(x, "fsvsim")) stop("Must be used on 'fsvsim' objects.")

 if (!length(i) == 1 || !is.numeric(i) || i < 1 || i > ncol(x$y))
  stop("Argument 'i' must be a single integer between 1 and ncol(x$y).")
 
 if (!length(j) == 1 || !is.numeric(j) || j < 1 || j > ncol(x$y))
  stop("Argument 'j' must be a single integer between 1 and ncol(x$y).")
 
 if (!is.numeric(these) || min(these) < 1 || max(these) > nrow(x$y))
  stop("Illegal argument value 'these'.")

 covelement <- (rep(x$facload[i,], each = length(these)) * exp(x$facvol[these,])) %*% x$facload[j,]
 if (i == j) covelement <- covelement + exp(x$idivol[these,i])
 
 as.numeric(covelement)
}


#' Extract "true" model-implied correlations of two series only
#'
#' \code{corelement} extracts the model-implied (time-varying) correlations between
#' (exactly) two component series.
#'
#' @param x Object of class \code{'fsvsim'}, usually resulting from a call
#' of the function \code{\link{fsvsim}}.
#' @param i Index of component series 1.
#' @param j Index of component series 2.
#' @param these Vector indicating which points in time should be extracted.
#'
#' @return Vector with the requested correlations.
#'
#' @family simulation
#'
#' @export

corelement <- function(x, i, j, these = seq_len(nrow(x$y))) {
 covij <- covelement(x, i, j, these)
 covii <- covelement(x, i, i, these)
 covjj <- covelement(x, j, j, these)
 covij / sqrt(covii * covjj)
}
