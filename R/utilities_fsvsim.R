#  #####################################################################################
#  R package factorstochvol by
#     Gregor Kastner Copyright (C) 2016-2021
#     Darjus Hosszejni Copyright (C) 2019-2021
#     Luis Gruber Copyright (C) 2021
#
#  This file is part of the R package factorstochvol: Bayesian Estimation
#  of (Sparse) Latent Factor Stochastic Volatility Models
#
#  The R package factorstochvol is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 2 or any
#  later version of the License.
#
#  The R package factorstochvol is distributed in the hope that it will
#  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the R package factorstochvol. If that is not the case,
#  please refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################

#' Extract "true" model-implied covariance matrix for several points in time
#'
#' \code{covmat} extracts the model-implied (time-varying) covariance matrix
#' from an \code{fsvsim} object.
#'
#' @param x Object of class \code{'fsvsim'}, usually resulting from a call
#' of the function \code{\link{fsvsim}}.
#' @param timepoints Vector indicating at which point(s) in time the
#' correlation matrices should be extracted. Can also be "all" or "last".
#' @param ... Ignored.
#'
#' @note Currently crudely implemented as an R loop over all time points,
#' may be slow.
#'
#' @return Array of dimension \code{m} times \code{m} times
#' \code{length(timepoints)}, containing the model-implied covariance matrix.
#'
#' @family simulation
#'
#' @export

covmat.fsvsim <- function(x, timepoints = "all", ...) {
 if (!is(x, "fsvsim")) stop("Argument 'x' must be of class 'fsvsim'.")

 if (is.character(timepoints)) {
   if (timepoints == "all") timepoints <- seq_len(ncol(x$fac)) else if (timepoints == "last") timepoints <- length(ncol(x$fac))
 } else if (!is.numeric(timepoints) || max(timepoints) > ncol(x$fac) || min(timepoints) < 1L) {
   stop("Illegal value for 'timepoints'.")
 }

 m <- nrow(x$facload)
 r <- ncol(x$facload)
 covmat <- array(NA_real_, dim = c(m, m, length(timepoints)))
 facload <- x$facload

 for (j in seq_along(timepoints)) {
  facvar <- exp(x$facvol[timepoints[j],])
  idivar <- exp(x$idivol[timepoints[j],])
  covmat[,,j] <- tcrossprod(sweep(facload, 2, facvar, '*'), facload)
  diag(covmat[,,j]) <- diag(covmat[,,j]) + idivar
 }

 covmat
}


#' Extract "true" model-implied correlation matrix for several points in time
#'
#' \code{cormat} extracts the model-implied (time-varying) covariance matrix
#' from an \code{fsvsim} object.
#'
#' @param x Object of class \code{'fsvsim'}, usually resulting from a call
#' of the function \code{\link{fsvsim}}.
#' @param timepoints Vector indicating at which point(s) in time the
#' correlation matrices should be extracted. Can also be "all" or "last".
#' @param ... Ignored.
#'
#' @note Currently crudely implemented as an R loop over all time points,
#' may be slow.
#'
#' @return Array of dimension \code{m} times \code{m} times
#' \code{length(timepoints)}, containing the model-implied correlation matrix.
#'
#' @family simulation
#'
#' @export

cormat.fsvsim <- function(x, timepoints = "all", ...) {
 mycovmat <- covmat.fsvsim(x, timepoints, ...)
 array(apply(mycovmat, 3, cov2cor), dim = dim(mycovmat))
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
