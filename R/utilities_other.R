#  #####################################################################################
#  R package factorstochvol by
#     Gregor Kastner Copyright (C) 2016-2020
#     Darjus Hosszejni Copyright (C) 2019-2020
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

#' Computes the log returns of a vector-valued time series
#'
#' \code{logret} computes the log returns of a multivariate time
#' series, with optional de-meaning.
#'
#' @param dat The raw data, a matrix or data frame with \code{n}
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
#' @rdname logret
#' @name logret
#' @export
logret.matrix <- function(dat, demean = FALSE, standardize = FALSE, ...) {
 tmp <- dat[,colSums(is.na(dat)) <= 0.5]
 tmp <- diff(log(as.matrix(tmp)))
 if (all(isTRUE(demean))) tmp <- tmp - rep(colMeans(tmp), each = nrow(tmp))
 if (all(isTRUE(standardize))) tmp <- tmp / rep(apply(tmp, 2, sd), each = nrow(tmp))
 tmp
}


#' @rdname logret
#' @name logret
#' @export
logret.data.frame <- function(dat, demean = FALSE, standardize = FALSE, ...) {
 dat <- data.matrix(dat)
 logret(dat, demean, standardize, ...)
}

#' @rdname logret
#' @name logret
#' @export
logret <- function(dat, demean = FALSE, standardize = FALSE, ...) {
 UseMethod("logret")
}



#' Ledermann bound for the number of factors
#' 
#' In the static factor case, the Ledermann bound is the largest
#' integer rank for which a unique decomposition of the covariance
#' matrix is possible. (This is the largest possible number of
#' factors which can be used for \code{\link[stats]{factanal}}.
#' 
#' @param m Number of component series.
#'
#' @return The Ledermann bound, a nonnegative integer.
#'
#' @seealso preorder
#'
#' @export
ledermann <- function(m) {
 as.integer(floor((2*m+1)/2 - sqrt((2*m+1)^2/4 - m^2 + m)))
}


#' Ad-hoc methods for determining the order of variables
#'
#' In factor SV models, the ordering of variables is often
#' chosen through a preliminary static factor analysis. These
#' methods are implemented in \code{preorder}.
#' After a maximum likelihood factor model fit to the data,
#' factor loadings are ordered as follows: The variable with the
#' highest loading on factor 1 is placed first, the variable with
#' the highest loading on factor 2 second (unless this variable
#' is already placed first, in which case the variable with the
#' second highest loading is taken).
#' 
#' @param dat Matrix containing the data, with \code{n} rows
#' (points in time) and \code{m} columns (component series).
#' @param factors Number of factors to be used, defaults to the
#' Ledermann bound.
#' @param type Can be "fixed" or "dynamic". The option "fixed"
#' means that that a \code{factors}-factor model is fit once and
#' the entire ordering is determined according to this fit
#' (the default). The option "dynamic" means that 
#' the model is re-fit \code{factors} times with the number of
#' factors going from 1 to
#' \code{factors} and in each round the correspondingly largest
#' loading is chosen.
#' @param transload Function for transforming the estimated
#' factor loadings before ordering. Defaults to the identity
#' function.
#'
#' @return A vector of length \code{m} with the ordering found.
#'
#' @seealso ledermann 
#'
#' @export
preorder <- function(dat, factors = ledermann(ncol(dat)), type = "fixed", transload = identity) {
 m <- ncol(dat)
 control <- list(opt = list(maxit = 100000)) 
 ordering <- rep(NA_integer_, m)
 
 if (type == "fixed") {
  fa <- factanal(dat, factors, control = control)
  for (i in 1:factors) {
   tmp <- order(transload(fa$loadings[,i]), decreasing = TRUE)
   ordering[i] <- tmp[!(tmp %in% ordering)][1]
  }
 } else if (type == "dynamic") {
  for (i in 1:factors) {
   fa <- factanal(dat, i, control = control)
   tmp <- order(transload(fa$loadings[,i]), decreasing = TRUE)
   ordering[i] <- tmp[!(tmp %in% ordering)][1]
  }
 } else stop("Unknown type")
 ordering[(factors + 1) : m] <- (1:m)[!((1:m) %in% ordering)]
 ordering
}


#' Ad-hoc method for (weakly) identifying the factor
#' loadings matrix
#'
#' In factor SV models, the identification of the factor loadings
#' matrix is often
#' chosen through a preliminary static factor analysis. 
#' After a maximum likelihood factor model is fit to the data,
#' variables are ordered as follows: The variable with the
#' lowest loadings on all factors except the first (relative to
#' it) is determined to lead the first factor,
#' the variable with the lowest loadings on all factors except the
#' first two (relative to these) is determined to lead the second
#' factor, etc. 
#' 
#' @param dat Matrix containing the data, with \code{n} rows
#' (points in time) and \code{m} columns (component series).
#' @param factors Number of factors to be used.
#' @param transload Function for transforming the estimated
#' factor loadings before ordering. Defaults to the absolute value
#' function.
#' @param relto Can be 'none', 'current' or 'all'. If 'none', the series
#' with the highest loadings is placed first, the series with the second
#' highest is placed second, and so on.
#' If 'current', the current factor loading is used as a reference, if 'all',
#' all previous loadings are summed up to be the reference.
#'
#' @return A \code{m} times \code{factors} matrix indicating
#' the restrictions.
#'
#' @note This function is automatically invoked by fsvsample if
#' restrict is set to 'auto'.
#'
#' @seealso ledermann 
#'
#' @export
findrestrict <- function(dat, factors, transload = abs, relto = 'all') {
 m <- ncol(dat)
 control <- list(opt = list(maxit = 100000)) 
 ordering <- rep(NA_integer_, m)
 restrict <- matrix(FALSE, nrow = m, ncol = factors)

 fa <- factanal(dat, factors, control = control)
 
 for (i in seq_len(factors - 1)) {
  if (relto == 'current') {
   cur <- transload(fa$loadings[, i]) 
   tmp <- order(rowSums(transload(fa$loadings[, (i+1):factors, drop = FALSE])) / cur)
  } else if (relto == 'all') {
   cur <- rowSums(transload(fa$loadings[, 1:i, drop = FALSE]))
   tmp <- order(rowSums(transload(fa$loadings[, (i+1):factors, drop = FALSE])) / cur)
  } else if (relto == 'none') {
   tmp <- order(transload(fa$loadings[,i]), decreasing = TRUE)
  } else stop("Unknown value of 'relto'.")
  ordering[i] <- tmp[!(tmp %in% ordering)][1]
  restrict[ordering[i], (i+1):factors] <- TRUE
 }
 ordering[factors : m] <- (1:m)[!((1:m) %in% ordering)]
 restrict
}


#' Computes the empirical exponentially weighted covariance matrix
#'
#' A common way to get estimates for time-varying covariance matrices
#' is the compute the exponentially weighted empirical covariance matrix.
#'
#' @param dat Matrix containing the data, with \code{n} rows
#' (points in time) and \code{m} columns (component series).
#'
#' @param alpha Speed of decay.
#' 
#' @param hist How far to go back in time?
#'
#' @return A \code{m} times \code{m} covariance matrix estimate.
#'
#' @export
expweightcov <- function(dat, alpha = 4/126, hist = 180) {
 n <- nrow(dat)
 mycov <- tcrossprod(dat[n - hist + 1,])
 for (i in 2:hist) {
  mycov <- (1 - alpha) * mycov + alpha * tcrossprod(dat[n - hist + i,])
 }
 mycov
}
