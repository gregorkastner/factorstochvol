rvolonly <- function(para, n) {
 
 mu <- para[1]
 phi <- para[2]
 sigma <- para[3]
 
 h <- rep(as.numeric(NA), n+1)
 h[1] <- rnorm(1, mean=mu, sd=sigma/sqrt((1-phi^2)))
 nu <- rnorm(n)

 # simulate w/ simple loop
 for (i in 2:(n+1)) h[i] <- mu + phi*(h[i-1]-mu) + sigma*nu[i-1]

 h
}


#' Simulate data from a factor SV model
#'
#' \code{fsvsim} generates simulated data from a factor SV model.
#'
#' @param n Length of the series to be generated.
#' @param series Number of component series \code{m}.
#' @param factors Number of factors \code{r}.
#' @param facload Can either be a matrix of dimension \code{m} times \code{r}
#' or one of the keywords "dense" and "sparse". If "dense" is chosen,
#' a (rather) dense lower triangular factor loadings matrix is randomly
#' generated. If "sparse" is chosen, a (rather) sparse lower triangular
#' factor loadings matrix is randomly generated.
#' @param idipara \emph{Optional} matrix of idiosyncratic SV parameters
#' to be used for simulation. Must have exactly three columns containing
#' the values of \code{mu}, \code{phi} and \code{sigma} for each
#' of \code{m} series, respectively. If omitted, plausible values are
#' generated.
#' @param facpara \emph{Optional} matrix of idiosyncratic SV parameters
#' to be used for simulation. Must have exactly two columns containing
#' the values of \code{phi} and \code{sigma} for each of \code{r} factors,
#' respectively. If omitted, plausible values are generated.
#' @param heteroskedastic Logical vector of length \code{m+r}. When
#' \code{TRUE}, time-varying volatilities are generated; when
#' \code{FALSE}, constant volatilities (equal to \code{mu}) are generated.
#' @param df If not equal to Inf, the factors are misspecified (come from
#' a t distribution instead of a Gaussian). Only used for testing.
#' 
#' @return The value returned is a list object of class \code{fsvsim} holding
#'  \itemize{
#'  \item{y}{The simulated data, stored in a \code{n} times \code{m} matrix with
#'  colnames 'Sim1', 'Sim2', etc.}
#'  \item{f}{The simulated factors, stored in a \code{r} times \code{r} matrix.}
#'  \item{facload}{Factor loadings matrix.}
#'  \item{facvol}{Latent factor log-variances for times 1 to \code{n}.}
#'  \item{facvol0}{Initial factor log-variances for time 0.}
#'  \item{facpara}{The parameters of the factor volatility processes.}
#'  \item{idivol}{Latent idiosyncratic log-variances for times 1 to \code{n}.}
#'  \item{idivol0}{Initial idiosyncratic log-variances for time 0.}
#'  \item{idipara}{The parameters of the idiosyncratic volatility
#'  processes.}
#' }
#' 
#' @note This object can be passed to many plotting functions to indicate
#' the data generating processes when visualizing results.
#' 
#' @export

fsvsim <- function(n = 1000, series = 10, factors = 1, facload = "dense", idipara, facpara,
		   heteroskedastic = rep(TRUE, series + factors), df = Inf) {
 if (!is.numeric(factors) || is.na(factors) || factors < 0) stop('Number of factors must be numeric and >= 0')
 if (!is.numeric(series) || is.na(series) || series < factors) stop('Number of series must be numeric and >= factors')
 
 if (length(facload) == 1 && is.character(facload)) {
  if (facload == "dense") {
   facload <- matrix(NA_real_, nrow = series, ncol = factors)
   if (factors >= 1) facload[,1] <- c(1, seq(from = 0.9, to = 0.1, length.out = series - 1))
   if (factors >= 2) facload[,2] <- c(0, 1, seq(from = 0.1, to = 0.8, length.out = series - 2))
   if (factors >= 3) facload[,3] <- c(0, 0, 1, seq(from = 0.7, to = 0.4, length.out = series - 3))
   if (factors >= 4) for (i in 4:factors) facload[,i] <- c(rep(0, i-1), 1, runif(series - i, -0.2, .8))
  } else if (facload == "sparse") {
   cutoff <- .2
   facload <- matrix(0, nrow = series, ncol = factors)
   for (i in 1:factors) {
    while (sum(abs(facload[,i]) > cutoff) < 3) {  # make sure each column has at least 3 nonzero elements (this is probably not enough for ident, but something at least...)
     facload[i,i] <- runif(1, cutoff + .1, 2 - i/factors)
     facload[(i+1):series,i] <- runif(series - i, -2*cutoff + cutoff*(i/factors), 1 - (1-cutoff-.25)*(i/factors))
    }
   }
   facload[abs(facload) < cutoff] <- 0
  } else stop('Strings allowed for "facload" are "dense" and "sparse"')
 } else {
  if (!is.matrix(facload)) stop('Factor loadings matrix "facload" must be a matrix, or a character vector equal to "dense", or "sparse"')
 }

 if (missing(idipara)) {
  idipara <- matrix(NA_real_, nrow = nrow(facload), ncol = 3)
  colnames(idipara) <- c("mu", "phi", "sigma")
  idipara[, "mu"] <- seq(from = -2, to = -1.1, length.out = series)
  idipara[, "phi"] <- seq(from = 0.8, to = 0.98, length.out = series)
  idipara[, "sigma"] <- seq(from = 0.6, to = 0.15, length.out = series)
 }
 if (!is.matrix(idipara)) stop('SV-parameter specification for idiosyncratic variances "idipara" must be a matrix')
 if (ncol(idipara) != 3) stop("idipara needs exactly three columns: mu, phi, sigma")
 if (nrow(idipara) != nrow(facload)) stop("Dimensions of idipara and facload don't match")
 
 if (missing(facpara)) {
  facpara <- matrix(NA_real_, nrow = ncol(facload), ncol = 2)
  colnames(facpara) <- c("phi", "sigma")
  facpara[, "phi"] <- c(0.99, 0.95, 0.97, runif(max(0, factors - 3), .95, .99))[seq_len(factors)]
  facpara[, "sigma"] <- c(0.1, 0.3, 0.1, runif(max(0, factors - 3), 0.1, 0.3))[seq_len(factors)]
 }
 if (!is.matrix(facpara)) stop('SV-parameter specification for factor variances "facpara" must be a matrix')
 if (ncol(facpara) != 2) stop("facpara needs exactly two columns: phi, sigma")
 if (nrow(facpara) != ncol(facload)) stop("Dimensions of facpara and facload don't match")

 # Some error checking for heteroskedastic
 if (length(heteroskedastic) != series + factors) stop("Argument 'heteroskedastic' must be of length series + factors.")
 if (!is.logical(heteroskedastic)) stop("Argument 'heteroskedastic' must be a vector containing only logical values.")
 if (is.null(heteroskedastic)) heteroskedastic <- rep(TRUE, series + factors)

 # simulate idiosyncratic variances:
 idivol <- matrix(NA_real_, nrow = n + 1, ncol = series)
 for (i in 1:series) {
  if (isTRUE(heteroskedastic[i])) {
   idivol[,i] <- rvolonly(idipara[i,], n)
  } else {
   idivol[,i] <- idipara[i,"mu"]
  }
 }
 
 # simulate factor variances and data:
 facvol <- matrix(NA_real_, nrow = n + 1, ncol = factors)
 if (factors > 0) {
  for (i in 1:factors) {
   if (isTRUE(heteroskedastic[series + i])) {
    facvol[,i] <- rvolonly(c(0, facpara[i,]), n)
   } else {
    facvol[,i] <- 0
   }
  }
  
  facvol0 <- facvol[1,,drop=FALSE] 
  
  # note: if df != Inf, this means misspecification!
  f <- t(apply(facvol[-1,,drop=FALSE], 2, function(x) exp(x/2) * rt(n, df = df)))
  #f <- t(apply(facvol[-1,,drop=FALSE], 2, function(x) exp(x/2) * runif(n, -2, 2)))
  tmp <- t(apply(idivol[-1,,drop=FALSE], 2, function(x) rnorm(n, mean=0, sd=exp(x/2))))
  y <- t(facload%*%f + tmp)
  
 } else {
  facvol0 <- matrix(NA_real_, nrow = 1, ncol = 0)
  f <- t(facvol)
  y <- apply(idivol[-1,,drop=FALSE], 2, function(x) rnorm(n, mean=0, sd=exp(x/2)))
 }

 colnames(y) <- paste0("Sim", 1:nrow(facload))

 ret <- list(y = y, f = f, facload = facload, facvol = facvol[-1,,drop=FALSE],
             facvol0 = facvol0, facpara = facpara, idivol = idivol[-1,,drop=FALSE],
             idivol0 = idivol[1,,drop=FALSE], idipara = idipara)

 class(ret) <- "fsvsim"
 ret
}
