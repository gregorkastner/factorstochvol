#' Pretty printing of an fsvsdraws object
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' of \code{\link{fsvsample}}.
#' @param ... Ignored.
#' 
#' @return Returns \code{x} invisibly.
#'
#' @family printing
#'
#' @export

print.fsvdraws <- function(x, ...) {
 m <- nrow(x$facload)
 r <- ncol(x$facload)
 n <- nrow(x$y)
 cat(paste("\nFitted factor stochastic volatility object with\n",
	   " -", formatC(m, width = 7), "series\n",
	   " -", formatC(r, width = 7), "factor(s)\n",
	   " -", formatC(n, width = 7), "timepoints\n",
	   " -", formatC(x$config$draws, width = 7), "MCMC draws\n",
	   " -", formatC(x$config$thin, width = 7), "thinning\n",
	   " -", formatC(x$config$burnin, width = 7), "burn-in\n\n"))
 invisible(x)
}


#' Extract summary statistics for the posterior covariance matrix
#' which have been stored during sampling
#'
#' \code{runningcovmat} extracts summary statistics from the model-implied
#' covariance matrix
#' from an \code{fsvdraws} object for one point in time.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' of \code{\link{fsvsample}}.
#' @param i A single point in time.
#' @param statistic Indicates which statistic should be extracted. Defaults
#' to \code{'mean'}.
#' @param type Indicates whether covariance (\code{cov}) or correlation
#' (\code{cor}) should be extracted.
#' 
#' @return Matrix containing the requested covariance matrix summary statistic.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(n = 500, series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1, runningstore = 6) # estimate
#'
#' cov100mean <- runningcovmat(res, 100) # extract mean at t = 100
#' cov100sd <- runningcovmat(res, 100, statistic = "sd") # extract sd
#' lower <- cov100mean - 2*cov100sd
#' upper <- cov100mean + 2*cov100sd
#'
#' true <- covmat(sim, 100) # true value
#'
#' # Visualize mean +/- 2sd and data generating values
#' par(mfrow = c(3,3), mar = c(2, 2, 2, 2))
#' for (i in 1:3) {
#'  for (j in 1:3) {
#'   plot(cov100mean[i,j], ylim = range(lower, upper), pch = 3,
#'   main = paste(i, j, sep = ' vs. '), xlab = '', ylab = '')
#'   lines(c(1,1), c(lower[i,j], upper[i,j]))
#'   points(true[i,j,1], col = 3, cex = 2)
#'  }
#' }
#' }
#'
#' @family extractors
#' 
#' @export


runningcovmat <- function(x, i, statistic = "mean", type = "cov") {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
  if (!exists("runningstore", x) || !exists(type, x$runningstore))
   stop("What you are requesting (argument 'type') hasn't been stored during sampling.")

 m <- ncol(x$y)
 n <- nrow(x$y)

 if (!length(i) == 1 || !is.numeric(i) || i < 1 || i > n)
  stop("Argument 'i' must be a single integer between 1 and ncol(x$y).")
 
 tryCatch(myobj <- x$runningstore[[type]][i,,statistic], error = function(e)
	   stop(paste0("Argument 'statistic' must be one of: ",
            paste(dimnames(x$runningstore[[type]])[[3]], collapse = ', '), ".")))

 ret <- diag(m)
 if (statistic != 'mean') diag(ret) <- 0
 if (type == "cor") {
  ret[lower.tri(ret)] <- myobj
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
 } else if (type == "cov") {
  ret[lower.tri(ret, diag = TRUE)] <- myobj
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
 } else stop("Argument 'type' must be 'cor' or 'cov'.")
 ret
}


#' Extract summary statistics for the posterior correlation matrix
#' which have been stored during sampling
#'
#' \code{runningcormat} extracts summary statistics from the model-implied
#' correlation matrix
#' from an \code{fsvdraws} object for one point in time.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' of \code{\link{fsvsample}}.
#' @param i A single point in time.
#' @param statistic Indicates which statistic should be extracted. Defaults
#' to \code{'mean'}.
#' @param type Indicates whether covariance (\code{cov}) or correlation
#' (\code{cor}) should be extracted.
#' 
#' @return Matrix containing the requested correlation matrix summary statistic.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(n = 500, series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1, runningstore = 6) # estimate
#'
#' cor100mean <- runningcormat(res, 100) # extract mean at t = 100
#' cor100sd <- runningcormat(res, 100, statistic = "sd") # extract sd
#' lower <- cor100mean - 2*cor100sd
#' upper <- cor100mean + 2*cor100sd
#'
#' true <- cov2cor(covmat(sim, 100)[,,1]) # true value
#'
#' # Visualize mean +/- 2sd and data generating values
#' par(mfrow = c(3,3), mar = c(2, 2, 2, 2))
#' for (i in 1:3) {
#'  for (j in 1:3) {
#'   plot(cor100mean[i,j], ylim = range(lower, upper), pch = 3,
#'   main = paste(i, j, sep = ' vs. '), xlab = '', ylab = '')
#'   lines(c(1,1), c(lower[i,j], upper[i,j]))
#'   points(true[i,j], col = 3, cex = 2)
#'  }
#' }
#' }
#'
#' @family extractors
#' 
#' @export

runningcormat <- function(x, i, statistic = "mean", type = "cor") {
 runningcovmat(x = x, i = i, statistic = statistic, type = type)
}


#' Extract posterior draws of the model-implied covariance matrix
#'
#' \code{covmat} extracts draws from the model-implied covariance matrix
#' from an \code{fsvdraws} object for all points in time which have been
#' stored.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' of \code{\link{fsvsample}}.
#' @param ... Ignored.
#' 
#' @note Currently crudely implemented as a double loop in pure R,
#' may be slow.
#' 
#' @return Array of dimension \code{m} times \code{m} times \code{draws}
#' times \code{timepoints} containing the posterior draws for the
#' model-implied covariance matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(n = 500, series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1, keeptime = "all") # estimate
#' covs <- covmat(res) # extract
#'
#' # Trace plot of determinant of posterior covariance matrix
#' # at time t = n = 500:
#' detdraws <- apply(covs[,,,500], 3, det)
#' ts.plot(detdraws)
#' abline(h = mean(detdraws), col = 2)          # posterior mean
#' abline(h = median(detdraws), col = 4)        # posterior median
#' abline(h = det(covmat(sim)[,,500]), col = 3) # implied by DGP
#'
#' # Trace plot of draws from posterior covariance of Sim1 and Sim2 at
#' # time t = n = 500:
#' ts.plot(covs[1,2,,500])
#' abline(h = covmat(sim)[1,2,500], col = 3) # "true" value
#' 
#' # Smoothed kernel density estimate:
#' plot(density(covs[1,2,,500], adjust = 2))
#'
#' # Summary statistics:
#' summary(covs[1,2,,500])
#' }
#'
#' @family extractors
#'
#' @seealso covmat covmat.fsvsim
#' 
#' @export

covmat.fsvdraws <- function(x, ...) {
 if (!is(x, "fsvdraws")) stop("Argument 'x' must be of class 'fsvdraws'.")
 m <- nrow(x$facload)
 r <- ncol(x$facload)
 draws <- dim(x$facload)[3]
 tpoints <- nrow(x$h)
 covmat <- array(NA_real_, dim = c(m, m, draws, tpoints))
 for (j in 1:tpoints) {
  for (i in 1:draws) {
   facload <- matrix(x$facload[,,i], nrow = m)
   facvar <- exp(x$h[j, m+(1:r),i])
   idivar <- exp(x$h[j, 1:m,i])
   covmat[,,i,j] <- tcrossprod(sweep(facload, 2, facvar, '*'), facload)
   diag(covmat[,,i,j]) <- diag(covmat[,,i,j]) + idivar
  }
 }
 covmat
}


#' Predicts factor and idiosyncratic log-volatilities h
#'
#' \code{predh} simulates from the posterior predictive distribution
#' of the latent log-variances h, both for factors as well as for 
#' idiosyncratic series.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#'
#' @return List of class \code{fsvpredh} containing two elements:
#' \itemize{
#' \item{idih}{Array containing the draws of the latent idiosyncratic
#' log-volatilities.}
#' \item{factorh}{Array containing the draws of the latent factor
#' log-volatilities.}
#' }
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1) # estimate
#' 
#' # Predict 1, 10, and 100 days ahead:
#' predobj <- predh(res, ahead = c(1, 10, 100))
#'
#' # Trace plot of draws from posterior predictive factor log-variance
#' # (one, ten, and 100 days ahead):
#' plot.ts(predobj$factorh[1,,])
#' 
#' # Smoothed kernel density estimates of predicted volas:
#' plot(density(exp(predobj$factorh[1,,"1"]/2), adjust = 2))
#' lines(density(exp(predobj$factorh[1,,"10"]/2), adjust = 2), col = 2)
#' lines(density(exp(predobj$factorh[1,,"100"]/2), adjust = 2), col = 3)
#' }
#'
#' @family predictors
#' 
#' @export

predh <- function(x, ahead = 1, each = 1) {
 if (!is(x, "fsvdraws")) stop("Argument 'x' must be of class 'fsvdraws'.")
 if (!is.vector(ahead) | !is.numeric(ahead) | any(is.na(ahead)))
  stop("Argument 'ahead' must be a numeric vector, NAs are not allowed.")
 ahead <- as.integer(ahead)
 ahead <- sort(ahead)
 if (any(ahead < 1)) stop("All elements of 'ahead' must be greater or equal to 1.")
 each <- as.integer(each)
 if (length(each) > 1 || each < 1) stop("Argument 'each' must be greater or equal to 1.")

 m <- ncol(x$y)
 r <- nrow(x$f)

 mus <- matrix(rep(x$para["mu",,], each), nrow = m + r)
 phis <- matrix(rep(x$para["phi",,], each), nrow = m + r)
 sigmas <- matrix(rep(x$para["sigma",,], each), nrow = m + r)
 hpredtmp <- matrix(rep(x$h[nrow(x$h),,], each), nrow = m + r)

 len <- ncol(sigmas)
 hpreds <- array(NA_real_, dim = c(m+r, len, length(ahead)))
 dimnames(hpreds) <- list(NULL, NULL, ahead = ahead)
 
 storecounter <- 1L
 for (i in 1:max(ahead)) {
  hpredtmp <- mus + phis * (hpredtmp - mus) + sigmas * rnorm(len*(m+r))
  if (i %in% ahead) {
   hpreds[,,storecounter] <- hpredtmp
   storecounter <- storecounter + 1L
  }
 }
 if (r > 0) {
  ret <- list(idih = hpreds[1:m,,,drop=FALSE], factorh = hpreds[m+(1:r),,,drop=FALSE])
 } else {
  ret <- list(idih = hpreds[1:m,,,drop=FALSE], factorh = array(NA_real_, dim = c(0, dim(hpreds)[2], dim(hpreds)[3])))
 }
 class(ret) <- c("fsvpredh")
 ret
}

#' Predicts covariance matrix
#'
#' \code{predcov} simulates from the posterior predictive distribution
#' of the model-implied covariance matrix.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#'
#' @return 4-dimensional array containing draws from the predictive
#' covariance distribution.
#' 
#' @note Currently crudely implemented as a triple loop in pure R,
#' may be slow.
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1) # estimate
#' 
#' # Predict 1, 10, and 100 days ahead:
#' predobj <- predcov(res, ahead = c(1, 10, 100))
#'
#' # Trace plot of draws from posterior predictive distribution
#' # of the covariance of Sim1 and Sim2:
#' # (one, ten, and 100 days ahead):
#' plot.ts(predobj[1,2,,])
#' 
#' # Smoothed kernel density estimates of predicted covariance
#' # of Sim1 and Sim2:
#' plot(density(predobj[1,2,,"1"], adjust = 2))
#' lines(density(predobj[1,2,,"10"], adjust = 2), col = 2)
#' lines(density(predobj[1,2,,"100"], adjust = 2), col = 3)
#' }
#'
#' @family predictors
#' 
#' @export

predcov <- function(x, ahead = 1, each = 1) {
 pred <- predh(x, ahead, each)
 m <- ncol(x$y)
 r <- nrow(x$f)
 ret <- array(NA_real_, dim = c(m, m, dim(pred$idih)[2], length(ahead)))
 dimnames(ret) <- list(colnames(x$y), colnames(x$y), NULL, ahead)
 for (i in seq_along(ahead)) {
  for (j in seq_len(dim(x$facload)[3])) {
   for (k in seq_len(each)) {
    tmp <- (j-1)*each+k
    if (r >= 1) {
     ret[,,tmp,i] <- x$facload[,,j] %*%
                    (exp(pred$factorh[,tmp,i]) * t(x$facload[,,j])) + 
                    diag(exp(pred$idih[,tmp,i]))
    } else {
     ret[,,tmp,i] <- diag(exp(pred$idih[,tmp,i]))
    }
   }
  }
 }
 ret
}

#' Predicts correlation matrix
#'
#' \code{predcor} simulates from the posterior predictive distribution
#' of the model-implied correlation matrix.
#' 
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#'
#' @return 4-dimensional array containing draws from the predictive
#' correlation distribution.
#' 
#' @note Currently crudely implemented as a triple loop in pure R,
#' may be slow.
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(series = 3, factors = 1) # simulate 
#' res <- fsvsample(sim$y, factors = 1) # estimate
#' 
#' # Predict 1, 10, and 100 days ahead:
#' predobj <- predcor(res, ahead = c(1, 10, 100))
#'
#' # Trace plot of draws from posterior predictive distribution
#' # of the correlation of Sim1 and Sim2:
#' # (one, ten, and 100 days ahead):
#' plot.ts(predobj[1,2,,])
#' 
#' # Smoothed kernel density estimates of predicted covariance
#' # of Sim1 and Sim2:
#' plot(density(predobj[1,2,,"1"], adjust = 2))
#' lines(density(predobj[1,2,,"10"], adjust = 2), col = 2)
#' lines(density(predobj[1,2,,"100"], adjust = 2), col = 3)
#' }
#'
#' @family predictors
#' 
#' @export

predcor <- function(x, ahead = 1, each = 1) {
 pred <- predh(x, ahead, each)
 m <- ncol(x$y)
 r <- nrow(x$f)
 ret <- array(NA_real_, dim = c(m, m, dim(pred$idih)[2], length(ahead)))
 dimnames(ret) <- list(colnames(x$y), colnames(x$y), NULL, ahead)
 for (i in seq_along(ahead)) {
  for (j in seq_len(dim(x$facload)[3])) {
   for (k in seq_len(each)) {
    tmp <- (j-1)*each+k
    if (r >= 1) {
     ret[,,tmp,i] <- cov2cor(x$facload[,,j] %*%
                     (exp(pred$factorh[,tmp,i]) * t(x$facload[,,j])) + 
                     diag(exp(pred$idih[,tmp,i])))
    } else {
     ret[,,tmp,i] <- diag(exp(pred$idih[,tmp,i]))
    }
   }
  }
 }
 ret
}


#' Predicts precision matrix and its determinant (Woodbury variant)
#'
#' \code{predprecWB} simulates from the posterior predictive distribution
#' of the model-implied precision matrix and its determinant
#' using the Woodbury matrix identity and the matrix determinant lemma
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#'
#' @return List containing two elements:
#' \itemize{
#' \item{precision}{Array containing the draws of the predicted
#' precision matrix.}
#' \item{precisionlogdet}{Matrix containing the draws of the determinant
#' of the predicted precision matrix.}
#' }
#' 
#' @note Currently crudely implemented as a triple loop in pure R,
#' may be slow.
#'
#' @family predictors
#' 
#' @seealso Usually used for evaluating the predictive likelihood when many
#' series but few factors are used, see
#' \code{\link{predloglik}} and \code{\link{predloglikWB}}.
#'
#' @export

predprecWB <- function(x, ahead = 1, each = 1) {
 pred <- predh(x, ahead, each)
 m <- ncol(x$y)
 r <- nrow(x$f)
 ret <- array(NA_real_, dim = c(m, m, dim(pred$idih)[2], length(ahead)))
 logdet <- array(NA_real_, dim = c(dim(pred$idih)[2], length(ahead)))
 dimnames(ret) <- list(colnames(x$y), colnames(x$y), NULL, ahead)
 dimnames(logdet) <- list(NULL, ahead)
 if (r > 0L) {
  for (i in seq_along(ahead)) {
   for (j in seq_len(dim(x$facload)[3])) {
    for (k in seq_len(each)) {
     tmp <- (j-1)*each+k
     idivars <- exp(pred$idih[,tmp,i])
     lala <- matrix(x$facload[,,j,drop=FALSE], ncol = r) / idivars
     lala2 <- diag(1/exp(pred$factorh[,tmp,i]), nrow = r) +
      crossprod(matrix(x$facload[,,j,drop=FALSE], ncol = r), lala)
     ret[,,tmp,i] <- diag(1/idivars, nrow = m) -
      lala %*% tcrossprod(solve(lala2), lala)
     logdet[tmp,i] <- determinant(lala2, logarithm = TRUE)$modulus +
      sum(pred$idih[,tmp,i]) + sum(pred$factorh[,tmp,i])
    }
   }
  }
 } else {
  for (i in seq_along(ahead)) {
   for (j in seq_len(dim(x$facload)[3])) {
    for (k in seq_len(each)) {
     tmp <- (j-1)*each+k
     ret[,,tmp,i] <- diag(1/exp(pred$idih[,tmp,i]), nrow = m) 
     logdet[tmp,i] <- sum(pred$idih[,tmp,i])
    }
   }
  }
 }
 list(precision = ret, precisionlogdet = 1/logdet)
}


#' Evaluates the predictive log likelihood using the predicted
#' covariance matrix
#'
#' \code{predloglik} approximates the predictive log likelihood by
#' simulating from the predictive distribution of the covariance
#' matrix and evaluating the corresponding multivariate normal
#' distribution.
#'
#' @param y Matrix of dimension \code{length(ahead)} times \code{m} where the
#' predictive density should be evaluated.
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#' @param alldraws Should all the draws be returned or just the final results?
#' (Can be useful to assess convergence.)
#'
#' @return Vector of length \code{length(ahead)} with log predictive
#' likelihoods.
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' # Simulate a time series of length 1100:
#' sim <- fsvsim(n = 1100, series = 3, factors = 1)
#' y <- sim$y
#'
#' # Estimate using only 1000 days:
#' res <- fsvsample(y[seq_len(1000),], factors = 1)
#' 
#' # Evaluate the 1, 10, and 100 days ahead predictive log
#' # likelihood:
#' ahead <- c(1, 10, 100)
#' scores <- predloglik(res, y[1000+ahead,], ahead = ahead, each = 10)
#' print(scores)
#' }
#'
#' @family predictors
#' 
#' @seealso Uses \code{\link{predcov}}. If \code{m} is large
#' but only few factors are used, consider also using
#' \code{\link{predloglikWB}}.
#'
#' @export

predloglik <- function(x, y, ahead = 1, each = 1, alldraws = FALSE) {
 predobj <- predcov(x, ahead, each)
 m <- ncol(x$y)
 r <- nrow(x$f)
 n <- dim(predobj)[3]
 if (!is.numeric(y) || !is.matrix(y) || ncol(y) != m || nrow(y) != length(ahead))
  stop("Argument 'y' must be a matrix of dimension c(length(ahead), nrow(x$y)).")
 ret <- array(NA_real_, dim = c(n, length(ahead)))
 realret <- rep(NA_real_, length(ahead))
 names(realret) <- ahead
 zeros <- matrix(0, nrow = m, ncol = n)
 for (i in seq_along(ahead)) {
  tmpy <- matrix(rep(y[i,], n), ncol = n)
  ret[,i] <- vecdmvnorm(tmpy, zeros, predobj[,,,i], log = TRUE)
  numericnormalizer <- max(ret[,i]) - 700 # exp(700) should be fine as double
  realret[i] <- log(mean(exp(ret[,i] - numericnormalizer))) + numericnormalizer
 }
 if (all(isTRUE(alldraws))) {
  return(list(predloglik = realret, predloglikdraws = ret))
 } else {
  return(realret)
 }
}


#' Evaluates the predictive log likelihood using the Woodbury identity
#'
#' \code{predloglikWB} approximates the predictive log likelihood exploiting
#' the factor structure and using the Woodbury idenitity and the
#' corresponding matrix determinant lemma. This is recommended only
#' if many series and few factors are present.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param y Matrix of dimension \code{length(ahead)} times \code{m} where the
#' predictive density should be evaluated.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#' @param alldraws Should all the draws be returned or just the final results?
#' (Can be useful to assess convergence.)
#'
#' @return Vector of length \code{length(ahead)} with log predictive
#' likelihoods.
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' # Simulate a time series of length 1100:
#' sim <- fsvsim(n = 1100, series = 3, factors = 1)
#' y <- sim$y
#'
#' # Estimate using only 1000 days:
#' res <- fsvsample(y[seq_len(1000),], factors = 1)
#' 
#' # Evaluate the 1, 10, and 100 days ahead predictive log
#' # likelihood:
#' ahead <- c(1, 10, 100)
#' scores <- predloglikWB(res, y[1000+ahead,], ahead = ahead, each = 10)
#' print(scores)
#' }
#'
#' @note Currently crudely implemented as a triple loop in pure R,
#' may be slow.
#'
#' @family predictors
#' 
#' @seealso Uses \code{\link{predprecWB}}. If \code{m} is small
#' or many factors are used, consider also using
#' \code{\link{predcov}}.
#'
#' @export

predloglikWB <- function(x, y, ahead = 1, each = 1, alldraws = FALSE) {
 covinvdet <- predprecWB(x, ahead, each)
 m <- ncol(x$y)
 r <- nrow(x$f)
 if (!is.numeric(y) || !is.matrix(y) || ncol(y) != m || nrow(y) != length(ahead))
  stop("Argument 'y' must be a matrix of dimension c(length(ahead), nrow(x$y)).")
 ret <- array(NA_real_, dim = c(dim(covinvdet$precisionlogdet)[1], length(ahead)))
 realret <- rep(NA_real_, length(ahead))
 names(realret) <- ahead
 for (i in seq_along(ahead)) {
  for (j in seq_len(dim(x$facload)[3])) {
   for (k in seq_len(each)) {
    tmp <- (j-1)*each+k
    ret[tmp,i] <- 1/covinvdet$precisionlogdet[tmp,i] +
                  crossprod(y[i,], covinvdet$precision[,,tmp,i]) %*% y[i,]
   }
  }
  ret[,i] <- -(m/2) * log(2*pi) - .5 * ret[,i]
  numericnormalizer <- max(ret[,i]) - 700 # exp(700) should be fine as double
  realret[i] <- log(mean(exp(ret[,i] - numericnormalizer))) + numericnormalizer
 }
 if (all(isTRUE(alldraws))) {
  return(list(predloglik = realret, predloglikdraws = ret))
 } else {
  return(realret)
 }
}


# OBSOLETE; is now implemented in C

predcondOLD <- function(x, ahead = 1, each = 1, ...) {
 if (!is(x, "fsvdraws")) stop("Argument 'x' must be of class 'fsvdraws'.")
 if (!is.vector(ahead) | !is.numeric(ahead) | any(is.na(ahead)))
  stop("Argument 'ahead' must be a numeric vector, NAs are not allowed.")
 ahead <- as.integer(ahead)
 ahead <- sort(ahead)
 if (any(ahead < 1)) stop("All elements of 'ahead' must be greater or equal to 1.")
 
 each <- as.integer(each)
 if (length(each) > 1 || each < 1) stop("Argument 'each' must be greater or equal to 1.")

 m <- ncol(x$y)
 r <- nrow(x$f)

 mus <- matrix(rep(x$para["mu",,], each), nrow = m + r)
 phis <- matrix(rep(x$para["phi",,], each), nrow = m + r)
 sigmas <- matrix(rep(x$para["sigma",,], each), nrow = m + r)
 hpredtmp <- matrix(rep(x$h[nrow(x$h),,], each), nrow = m + r)

 len <- ncol(sigmas)
 hpreds <- array(NA_real_, dim = c(m+r, len, length(ahead)))
 facpreds <- array(NA_real_, dim = c(r, len, length(ahead)))
 meanpreds <- array(NA_real_, dim = c(m, len, length(ahead)))
 
 dimnames(hpreds) <- dimnames(meanpreds) <- list(NULL, NULL, ahead = ahead)

 storecounter <- 1L
 for (i in 1:max(ahead)) {
  hpredtmp <- mus + phis * (hpredtmp - mus) + sigmas * rnorm(len*(m+r))
  if (i %in% ahead) {
   hpreds[,,storecounter] <- hpredtmp
   facpreds[,,storecounter] <- exp(hpredtmp[m+seq_len(r),]/2) * rnorm(len*r)
   if (r >= 1) {
    for (j in 1:len) {  # SLOW!
     meanpreds[,j,storecounter] <- matrix(x$facload[,,ceiling(j/each)], ncol=r) %*%
                                   facpreds[,j,storecounter] 
    }
   } else {
    meanpreds[,,storecounter] <- 0
   }
   storecounter <- storecounter + 1L
  }
 }

 ret <- list(meanpreds, exp(hpreds[1:m,,,drop=FALSE]/2))
 names(ret) <- c("means", "vols")
 class(ret) <- c("fsvpredcond")
 ret
}


#' A posteriori sign identification
#'
#' \code{signident} provides methods for identifying the signs of
#' the factor loadings after running the MCMC sampler
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param method Can be "diagonal" or "maximin". If "diagonal" is
#' chosen, the diagonal elements of the factor loadings matrix
#' are assumed to have positive signs
#' and the others are arranged accordingly.
#' If "maximin" is chosen, for each factor, \code{signident} looks
#' for the series where
#' the minimum absolute loadings are biggest and chooses this series
#' to have positive loadings.
#' @param implementation Either 1, 2, or 3 (the default). Determines
#' how the reordering is implemented. Should not be necessary to depart
#' from the default.
#' 
#' @return Returns an object of class \code{'fsvdraws'} with adjusted
#' factors and factor loadings. Moreover, a list element called
#' \code{'identifier'} is added, providing the numbers of the series
#' used for identification and the corresponding minimum distances to
#' zero.
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- fsvsim(series = 8, factors = 2) # simulate 
#' res <- fsvsample(sim$y, factors = 2, signswitch = TRUE,
#'                  draws = 2000, burnin = 1000) # estimate
#' 
#' # Plot unidentified loadings:
#' facloaddensplot(res, fsvsimobj = sim, rows = 8)
#'
#' # Identify:
#' res <- signident(res)
#'
#' # Plot identified loadings:
#' facloaddensplot(res, fsvsimobj = sim, rows = 8)
#' }
#'
#' @family postprocessing
#' 
#' @export

signident <- function(x, method = "maximin", implementation = 3) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (method != "diagonal" & method != "maximin")
  stop("Argument 'method' must either be 'diagonal' or 'maximin'.")
 r <- dim(x$facload)[2]
 ftpoints <- dim(x$f)[2]

 if (r == 0) {
  x$identifier <- matrix(NA_real_, nrow = 0, ncol = 2)
 } else {
  identifier <- 1:r  # if method == "diagonal", overwritten otherwise
  distance <- rep(NA_real_, r)

  for (i in 1:r) {
   faccol <- matrix(x$facload[,i,,drop=FALSE], nrow = nrow(x$facload))
   if (method == "maximin") { # for each factor, look for the series where the
                              # minimum absolute loadings are biggest
    identifier[i] <- which.max(apply(abs(faccol), 1, min))
   }

   distance[i] <- max(apply(abs(faccol), 1, min))
   mysig <- sign(faccol[identifier[i],])
   x$facload[,i,] <- t(t(faccol) * mysig)
   if (implementation == 1) {
    x$f[i,,] <- x$f[i,,] * rep(mysig, each = ftpoints)
   } else if (implementation == 2) {
    for (j in seq(along = mysig)) {
     x$f[i,,j] <- x$f[i,,j] * mysig[j]
    }
   } else if (implementation == 3) {
    for (j in 1:ftpoints) {
     x$f[i,j,] <- x$f[i,j,] * mysig
    }
   } else stop("Err0r.")
  }
  x$identifier <- matrix(c(identifier, distance), ncol = 2)
 }

 colnames(x$identifier) <- c("identifier", "distance")
 
 x
}


#' A posteriori factor order identification
#'
#' \code{orderident} provides some (very ad-hoc) methods for identifying
#' the ordering of the factors after running the (unrestricted) MCMC
#' sampler by 
#' ordering according to the argument \code{method}.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param method Methods currently supported:
#' \itemize{
#' \item \code{summean} Sort by sum of mean loadings (descending).
#' \item \code{summeaninv} Sort by sum of mean loadings (ascending).
#' \item \code{summeanabs} Sort by sum of mean absolute loadings (descending).
#' \item \code{summed} Sort by sum of median loadings (descending).
#' \item \code{summedinv} Sort by sum of median loadings (ascending).
#' \item \code{summedabs} Sort by sum of median absolute loadings (descending).
#' \item \code{maxmed} Sort by maximum median loadings (descending).
#' \item \code{maxmedinv} Sort by maximum median loadings (ascending).
#' \item \code{maxmedrel} Sort by maximum median loadings, relative to the sum of all median loadings on that factor (descending).
#' \item \code{maxmedabsrel} Sort by maximum absolute median loadings, relative to the sum of all median loadings on that factor (descending).
#' }
#'
#' @return Returns an object of class \code{'fsvdraws'} with adjusted
#' ordering.
#'
#' @family postprocessing
#' 
#' @export

orderident <- function(x, method = "summed") {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 possiblemethods <- c("summean", "summeaninv", "summeanabs", "summedabs", "summed",
		      "summedinv", "maxmed", "maxmedinv", "maxmedrel", "maxmedrelabs")
 if (is.numeric(method)) {
  if (length(method) != ncol(x$facload))
   stop(paste("Argument 'method' must be numeric of length 'number of factors' or one of:",
	      paste(possiblemethods, collapse = ", ")))
  myorder <- NULL
  for (i in seq_len(ncol(x$facload))) {
   tmp <- rev(order(apply(x$facload[method[i],,], 1, median)))
   for (j in seq_len(ncol(x$facload))) {
    if (!(tmp[j] %in% myorder)) {
     myorder[i] <- tmp[j]
     break
    }
   }
  }
 } else {
  if (!all(method %in% possiblemethods))
   stop(paste("Argument 'method' must be numeric of length 'number of factors' or one of:",
	      paste(possiblemethods, collapse = ", ")))
  r <- dim(x$facload)[2]
  if (r <= 1) return(x)
  m <- dim(x$facload)[1]
  orderstats <- switch(method,
    summean = colSums(apply(x$facload, 1:2, mean)),
    summeaninv = colSums(apply(x$facload, 1:2, mean)),
    summeanabs = colSums(apply(abs(x$facload), 1:2, mean)),
    summed = colSums(apply(x$facload, 1:2, median)),
    summedinv = colSums(apply(x$facload, 1:2, median)),
    summedabs = colSums(apply(abs(x$facload), 1:2, median)),
    maxmed = apply(apply(x$facload, 1:2, median), 2, max),
    maxmedinv = apply(apply(x$facload, 1:2, median), 2, max),
    maxmedrel = apply(apply(x$facload, 1:2, median), 2, max) / colSums(apply(x$facload, 1:2, median)),
    maxmedabsrel = apply(apply(abs(x$facload), 1:2, median), 2, max) / colSums(apply(abs(x$facload), 1:2, median)),
  )
  myorder <- order(orderstats, decreasing = TRUE)
  if (method %in% possiblemethods[grep("inv", possiblemethods)]) myorder <- rev(myorder)
 }
 m <- nrow(x$facload)
 r <- ncol(x$facload)
 x$facload <- x$facload[,myorder,,drop=FALSE]
 x$f <- x$f[myorder,,,drop=FALSE]
 x$para[,m+(1:r),] <- x$para[,m+myorder,,drop=FALSE]
 x$h[,m+(1:r),] <- x$h[,m+myorder,,drop=FALSE]
 if (exists("h", x$runningstore)) x$runningstore$h[,m+(1:r),] <- x$runningstore$h[,m+myorder,,drop=FALSE]
 if (exists("f", x$runningstore)) x$runningstore$f <- x$runningstore$f[myorder,,,drop=FALSE]
 if (exists("sd", x$runningstore)) x$runningstore$sd[,m+(1:r),] <- x$runningstore$sd[,m+myorder,,drop=FALSE]
 if (exists("identifier", x)) x$identifier <- x$identifier[myorder,]
 x
}
