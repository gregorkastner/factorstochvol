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

#' Markov Chain Monte Carlo (MCMC) Sampling for the Factor Stochastic
#' Volatility Model.
#'
#' \code{fsvsample} simulates from the joint posterior distribution and returns
#' the MCMC draws. It is the main workhorse to conduct inference for factor
#' stochastic volatility models in this package.
#'
#' @param y Data matrix. Each of \code{m} columns is assumed to contain a single
#' (univariate) series of length \code{n}.
#'
#' @param factors Number of latent factors to be estimated.
#'
#' @param draws Number of MCMC draws kept after burn-in.
#'
#' @param burnin Number of initial MCMC draws to be discarded.
#'
#' @param zeromean Logical. If \code{FALSE}, a constant mean is included in
#' the model for each of the \code{m} univariate series.
#' If \code{TRUE}, the mean is not modeled. Defaults to \code{TRUE}.
#'
#' @param priormu Vector of length 2 denoting prior mean and standard deviation
#' for unconditional levels of the idiosyncratic log variance processes.
#'
#' @param priorphiidi Vector of length 2, indicating the shape parameters for the
#' Beta prior distributions of the transformed parameters \code{(phi+1)/2}, where
#' \code{phi} denotes the persistence of the idiosyncratic log variances.
#'
#' @param priorphifac Vector of length 2, indicating the shape parameters for the
#' Beta prior distributions of the transformed parameters \code{(phi+1)/2}, where
#' \code{phi} denotes the persistence of the factor log variances.
#'
#' @param priorsigmaidi Vector of length \code{m} containing the
#' prior volatilities of log variances. If \code{priorsigmaidi} has exactly
#' one element, it will be recycled for all idiosyncratic log variances.
#'
#' @param priorsigmafac Vector of length \code{factors} containing the
#' prior volatilities of log variances. If \code{priorsigmafac} has exactly
#' one element, it will be recycled for all factor log variances.
#'
#' @param priorfacloadtype Can be \code{"normal"}, \code{"rowwiseng"},
#' \code{"colwiseng"}.
#' \itemize{
#'  \item{\code{"normal"}: }{Normal prior. The value of \code{priorfacload}
#'                           is interpreted as the standard deviations of the
#'                           Gaussian prior distributions for the factor loadings.}
#'  \item{\code{"rowwiseng"}: }{Row-wise Normal-Gamma prior. The value of \code{priorfacload}
#'                              is interpreted as the shrinkage parameter \code{a}.}
#'  \item{\code{"colwiseng"}: }{Column-wise Normal-Gamma prior. The value of \code{priorfacload}
#'                              is interpreted as the shrinkage parameter \code{a}.}
#' }
#' For details please see Kastner (2019).
#'
#' @param priorfacload Either a matrix of dimensions \code{m} times \code{factors}
#' with positive elements or a single number (which will be recycled accordingly).
#' The meaning of \code{priorfacload} depends on the setting of \code{priorfacloadtype}
#' and is explained there.
#'
#' @param facloadtol Minimum number that the absolute value of a factor loadings draw can take.
#' Prevents numerical issues that can appear when strong shrinkage is enforced if chosen
#' to be greater than zero.
#'
#' @param priorng Two-element vector with positive entries indicating the Normal-Gamma
#' prior's hyperhyperparameters \code{c} and \code{d}.
#'
#' @param priorh0idi Vector of length 1 or \code{m}, containing
#' information about the Gaussian prior for the initial idiosyncratic
#' log variances.
#' If an element of \code{priorh0idi} is a nonnegative number, the conditional
#' prior of the corresponding initial log variance h0 is assumed to be Gaussian
#' with mean 0 and standard deviation \code{priorh0idi} times $sigma$.
#' If an element of
#' \code{priorh0idi} is the string 'stationary', the prior of the corresponding
#' initial log volatility is taken to be from the stationary distribution, i.e.
#' h0 is assumed to be Gaussian with mean 0 and variance $sigma^2/(1-phi^2)$.
#'
#' @param priorh0fac Vector of length 1 or \code{factors}, containing
#' information about the Gaussian prior for the initial factor
#' log variances.
#' If an element of \code{priorh0fac} is a nonnegative number, the conditional
#' prior of the corresponding initial log variance h0 is assumed to be Gaussian
#' with mean 0 and standard deviation \code{priorh0fac} times $sigma$.
#' If an element of
#' \code{priorh0fac} is the string 'stationary', the prior of the corresponding
#' initial log volatility is taken to be from the stationary distribution, i.e.
#' h0 is assumed to be Gaussian with mean 0 and variance $sigma^2/(1-phi^2)$.
#'
#' @param priorbeta numeric vector of length 2, indicating the mean and
#' standard deviation of the Gaussian prior for the regression parameters. The
#' default value is \code{c(0, 10000)}, which constitutes a very vague prior
#' for many common datasets. Not used if \code{zeromean} is \code{TRUE}.
#'
#' @param thin Single number greater or equal to 1, coercible to integer.
#' Every \code{thin}th MCMC draw is kept and returned. The default value
#' is 1, corresponding to no thinning of the draws, i.e. every draw is
#' stored.
#'
#' @param keeptime Either a number coercible to a positive integer, or a string
#' equal to "all" or "last". If a number different from 1 is provided, only every
#' \code{keeptime}th latent log-volatility is being monitored. If, e.g.,
#' \code{keeptime = 3}, draws for the latent log variances
#' \code{h_1,h_4,h_7,...} will
#' be kept. If \code{keeptime} is set to "all", this is equivalent to setting it
#' to 1. If
#' \code{keeptime} is set to "last" (the default), only draws for the very last
#' latent log variances h_n are kept.
#'
#' @param runningstore Because most machines these days do not have enough memory
#' to store all draws for all points in time, setting \code{runningstore} to an
#' integer greater than 0 will cause \code{fsvsample} to store the first
#' \code{runningstoremoments}
#' ergodic moments of certain variables of interest. More specifically, mean,
#' variance, skewness, etc. will be stored for certain variables
#' if \code{runningstore} is set to a value...
#' \itemize{
#'  \item{\code{>= 1}: }{Latent log variances \code{h_1,h_2,...,h_(n+r)}.}
#'  \item{\code{>= 2}: }{Latent factors \code{f_1,...,f_r}.}
#'  \item{\code{>= 3}: }{Latent volatilities \code{sqrt(exp(h_1,h_2,...,h_(n+r)))}.}
#'  \item{\code{>= 4}: }{Conditional covariance matrix and the square roots of its
#'                      diagonal elements.}
#'  \item{\code{>= 5}: }{Conditional correlation matrix.}
#'  \item{\code{>= 6}: }{Communalities, i.e. proportions of variances explained
#'                      through the common factors.}
#' }
#'
#' @param runningstorethin How often should the calculation of running moments be
#' conducted? Set to a value > 1 if you want to avoid time consuming calculations
#' at every MCMC iteration.
#'
#' @param runningstoremoments Selects how many running moments (up to 4) should be
#' calculated.
#'
#' @param quiet Logical value indicating whether the progress bar and other
#' informative output during sampling should be omitted. The default
#' value is \code{FALSE}, implying verbose output.
#'
#' @param restrict Either "upper", "none", or "auto", indicating whether
#' the factor loadings
#' matrix should be restricted to have zeros above the diagonal ("upper"),
#' whether all elements should be estimated from the data ("none"), or
#' whether the function \code{\link{findrestrict}} should be invoked for a
#' priori finding suitable zeros. Setting
#' \code{restrict} to "upper" or "auto" often stabilizes MCMC
#' estimation and can be important for identifying the factor loadings matrix,
#' however, it generally is a strong prior assumption. Setting
#' \code{restrict} to "none" is usually the preferred option if identification
#' of the factor loadings matrix is of less concern but covariance estimation
#' or prediction is the goal. Alternatively, \code{restrict} can be a
#' logical matrix
#' of dimension \code{c(m, r)} indicating which elements should be unrestricted
#' (where \code{restrict} is \code{FALSE}) or zero (where \code{restrict} is
#' \code{TRUE}).
#'
#' @param interweaving The following values for interweaving the factor loadings
#' are accepted:
#' \itemize{
#'  \item{0: }{No interweaving.}
#'  \item{1: }{Shallow interweaving through the diagonal entries.}
#'  \item{2: }{Deep interweaving through the diagonal entries.}
#'  \item{3: }{Shallow interweaving through the largest absolute entries in each column.}
#'  \item{4: }{Deep interweaving through the largest absolute entries in each column.}
#' }
#' For details please see Kastner et al. (2017).
#' A value of 4 is the highly recommended default.
#'
#' @param signswitch Set to \code{TRUE} to turn on a random sign switch of
#' factors and loadings. Note that the signs of each factor loadings matrix column
#' and the corresponding factor cannot be identified from the likelihood.
#'
#' @param heteroskedastic Vector of length 1, 2, or \code{m + factors},
#' containing logical values indicating whether time-varying
#' (\code{heteroskedastic = TRUE}) or constant (\code{heteroskedastic = FALSE})
#' variance should be estimated.
#' If \code{heteroskedastic} is of length 2 it will be recycled accordingly,
#' whereby the first element is used for all idiosyncratic variances and
#' the second element is used for all factor variances.
#'
#' @param priorhomoskedastic Only used if at least one element of
#' \code{heteroskedastic} is set to \code{FALSE}. In that case,
#' \code{priorhomoskedastic} must be a matrix with positive entries
#' and dimension c(m, 2). Values in column 1 will be interpreted as the
#' shape and values in column 2 will be interpreted as the rate parameter
#' of the corresponding inverse gamma prior distribution of the idisyncratic
#' variances.
#'
#' @param expert \emph{optional} named list of expert parameters for the
#' univariate SV models (will be transformed and passed to the \code{stochvol} package). For most
#' applications, the default values probably work best. Interested users are
#' referred to Kastner and Frühwirth-Schnatter (2014), the package vignette, and Kastner (2016). If
#' \code{expert} is provided, it may contain the following named elements:
#' \itemize{
#'  \item{\code{parameterization}: }{Character string equal to \code{"centered"},
#'        \code{"noncentered"}, \code{"GIS_C"}, or \code{"GIS_NC"}. Defaults to
#'        \code{"GIS_C"}.}
#'  \item{\code{mhcontrol}: }{Single numeric value controlling the proposal
#'        density of a Metropolis-Hastings (MH) update step when sampling
#'        \code{sigma}. If \code{mhcontrol} is smaller than 0, an independence
#'        proposal will be used, while values greater than zero control the
#'        stepsize of a log-random-walk proposal. Defaults to \code{-1}.}
#'  \item{\code{gammaprior}: }{Single logical value indicating whether a Gamma
#'        prior for \code{sigma^2} should be used. If set to \code{FALSE},
#'        an Inverse Gamma prior is employed. Defaults to \code{TRUE}.}
#'  \item{\code{truncnormal}: }{Single logical value indicating whether a
#'        truncated Gaussian distribution should be used as proposal for
#'        draws of \code{phi}. If set to \code{FALSE}, a regular Gaussian
#'        prior is employed and the draw is immediately discarded when values
#'        outside the unit ball happen to be drawn. Defaults to \code{FALSE}.}
#'  \item{\code{mhsteps}: }{Either \code{1}, \code{2}, or \code{3}. Indicates
#'        the number of blocks used for drawing from the posterior of the
#'        parameters. Defaults to \code{2}.}
#'  \item{\code{proposalvar4sigmaphi}: }{Single positive number indicating the
#'        conditional prior variance of \code{sigma*phi} in the ridge
#'        \emph{proposal} density for sampling \code{(mu, phi)}.
#'        Defaults to \code{10^8}.}
#'  \item{\code{proposalvar4sigmatheta}: }{Single positive number indicating
#'        the conditional prior variance of \code{sigma*theta} in the ridge
#'        \emph{proposal} density for sampling \code{(mu, phi)}.
#'        Defaults to \code{10^12}.}
#' }
#'
#' @param startpara \emph{optional} numeric matrix of dimension
#' \code{c(3, m + factors)}, containing the starting values
#' for the parameter draws. The first \code{m} columns must contain
#' parameters values corresponding to the idiosyncratic volatilities,
#' the subsequent \code{factor} columns must contain parameter values
#' corresponding to the factor volatilities. The first row of \code{startpara}
#' corresponds to \code{mu}, the level of the log variances (can be arbitrary
#' numerical values), the second row corresponds to \code{phi}, the persistence
#' parameters of the log variances (numeric values between \code{-1} and \code{1}),
#' and the third row corresponds to \code{sigma} (positive numeric values).
#'
#' @param startlogvar \emph{optional} numeric matrix of dimension
#' \code{c(n, m + factors)}, containing the starting values of the
#' latent log variances.
#' The first \code{m} rows correspond to the idiosyncratic log variances,
#' the subsequent \code{factor} rows correspond to the factor log variances.
#' Was previously called startlatent.
#'
#' @param startlatent \emph{Deprecated.} Please use \code{startlogvar} instead.
#'
#' @param startlogvar0 \emph{optional} numeric vector of length
#' \code{m + factors}, containing the starting values of the initial latent
#' log variances.
#' The first \code{m} elements correspond to the idiosyncratic log variances,
#' the subsequent \code{factor} elements correspond to the factor log variances.
#' Was previously called startlatent0.
#'
#' @param startlatent0 \emph{Deprecated.} Please use \code{startlogvar0} instead.
#'
#' @param startfacload \emph{optional} numeric matrix of dimension
#' \code{c(m, factors)}, containing the starting values of the
#' factor loadings.
#' In case of a single factor model, a numeric vector of length \code{n} is also accepted.
#'
#' @param startfac \emph{optional} numeric matrix of dimension
#' \code{c(factors, n)}, containing the starting values of the
#' latent factors.
#' In case of a single factor model, a numeric vector of length \code{n} is also accepted.
#'
#' @param startfacloadvar \emph{optional} numeric matrix of dimension
#' \code{c(m, factors)}, containing the starting values of the
#' factor loadings variances \eqn{\tau_{ij}^2}. Used only when the normal-gamma
#' prior is employed (priorfacloadtype != "normal") while ignored when static
#' loadings variances are used (priorfacloadtype == "normal").
#'
#' @param samplefac If set to \code{FALSE}, the factors are not sampled (but
#' remain at their starting values forever). This might be useful if one
#' wants to include observed factors instead of latent ones.
#'
#' @param signident If set to \code{FALSE}, no ex-post sign-identification is
#' performed. Defaults to \code{TRUE}.
#'
#' @details For details concerning the factor SV algorithm please see
#' Kastner et al. (2017), details about the univariate SV estimation
#' can be found in Kastner and Frühwirth-Schnatter (2014).
#'
#' @return The value returned is a list object of class \code{fsvdraws} holding
#'  \itemize{
#'  \item{\code{facload}: }{Array containing draws from the posterior distribution of the
#'                 factor loadings matrix.}
#'  \item{\code{fac}: }{Array containing factor draws from the posterior distribution.}
#'  \item{\code{logvar}: }{Array containing idiosyncratic and factor initial log variance draws.}
#'  \item{\code{logvar0}: }{Array containing idiosyncratic and factor log variance draws.}
#'  \item{\code{para}: }{Array containing parameter draws form the posterior distribution.}
#'  \item{\code{y}: }{Matrix containing the data supplied.}
#'  \item{\code{latestauxiliary}: }{List containing the latest draws of auxiliary quantities used for
#'                      sampling the factor loadings matrix.}
#'  \item{\code{runningstore}: }{List whose elements contain ergodic moments of certain
#'                      variables of interest. See argument
#'                      \code{runningstore} for details about what
#'                      is being stored here.}
#'  \item{\code{config}: }{List containing information on configuration parameters.}
#'  \item{\code{priors}: }{List containing prior hyperparameter values.}
#'  \item{\code{identifier}: }{Matrix containing the indices of the series used for
#'    ex-post sign-identification along with the corresponding minimum distances to zero.
#'    See \code{\link{signident}} for details.}
#' }
#' To display the output, use \code{print}, \code{plot}, and in particular specialized
#' extractors and printing functions.
#' The \code{print} method prints a high-level overview; specialized extractors such as
#' \code{\link{covmat}} or \code{\link{runningcovmat}} are also available.
#' The \code{plot} method invokes a simple covariance matrix plot; specialized plotting
#' functions are linked in the documentation of \code{\link{plot.fsvdraws}}.
#'
#' @references Kastner, G., Frühwirth-Schnatter, S., and Lopes, H.F. (2017).
#' Efficient Bayesian Inference for Multivariate Factor Stochastic Volatility Models.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26}(4), 905--917,
#' \doi{10.1080/10618600.2017.1322091}.
#'
#' @references Kastner, G. (2019).
#' Sparse Bayesian Time-Varying Covariance Estimation in Many Dimensions
#' \emph{Journal of Econometrics}, \bold{210}(1), 98--115,
#' \doi{10.1016/j.jeconom.2018.11.007}
#'
#' @references Kastner, G. (2016).
#' Dealing with stochastic volatility in time series using the R package
#' stochvol.
#' \emph{Journal of Statistical Software}, \bold{69}(5), 1--30,
#' \doi{10.18637/jss.v069.i05}.
#'
#' @references Kastner, G. and Frühwirth-Schnatter, S. (2014).
#' Ancillarity-Sufficiency Interweaving Strategy (ASIS) for Boosting MCMC
#' Estimation of Stochastic Volatility Models.
#' \emph{Computational Statistics & Data Analysis}, \bold{76}, 408--423,
#' \doi{10.1016/j.csda.2013.01.002}.
#'
#' @family wrappers
#'
#' @examples
#' \donttest{
#' # Load exchange rate data (ships with stochvol):
#' data(exrates, package = "stochvol")
#' exrates$date <- NULL
#'
#' # Compute the percentage log returns:
#' dat <- 100 * logret(exrates)
#'
#' # We are going to fit a one-factor model so the ordering is irrelevant
#' # NOTE that these are very few draws, you probably want more...
#' res <- fsvsample(dat, factors = 2, draws = 2000, burnin = 1000,
#'   runningstore = 6, zeromean = FALSE)
#'
#' voltimeplot(res)
#'
#' corimageplot(res, nrow(dat), plotCI = 'circle')
#'
#' oldpar <- par(ask = TRUE)
#' plot(res)
#' par(oldpar)
#' pairs(t(res$beta[1:4, ]))
#' }
#'
#' @export
fsvsample <- function(y, factors = 1, draws = 1000, thin = 1, burnin = 1000,
                      restrict = "none", zeromean = TRUE,
                      priorfacloadtype = "rowwiseng", priorfacload = .1,
                      facloadtol = 1e-18,
                      priorng = c(1, 1), priormu = c(0, 10),
                      priorphiidi = c(10, 3), priorphifac = c(10, 3),
                      priorsigmaidi = 1, priorsigmafac = 1,
                      priorh0idi = "stationary", priorh0fac = "stationary",
                      priorbeta = c(0, 10000),
                      keeptime = "last", heteroskedastic = TRUE,
                      priorhomoskedastic = NA, runningstore = 6,
                      runningstorethin = 10, runningstoremoments = 2,
                      signident = TRUE, signswitch = FALSE, interweaving = 4,
                      quiet = FALSE, samplefac = TRUE, startfac, startpara,
                      startlogvar, startlatent, startlogvar0, startlatent0,
                      startfacload, startfacloadvar, expert) {

 # startlatent and startlogvar0 are being faded out
 if (!missing("startlatent") || !missing("startlatent0")) {
  warning("The arguments 'startlatent' and 'startlatent0' were renamed to 'startlogvar' and 'startlogvar0', respectively. Please use the new argument names.")
  if (missing("startlogvar") && !missing("startlatent")) startlogvar <- startlatent
  if (missing("startlogvar0") && !missing("startlatent0")) startlogvar0 <- startlatent0
 }

 # Some error checking for y
 if (is(y, "fsvsim")) {
  y <- y[["y"]]
  warning("Extracted data vector from 'fsvsim'-object.")
 }
 if (!is.matrix(y)) stop("Argument 'y' must be a matrix.")

 n <- nrow(y)
 m <- ncol(y)

 if (n < m) warning(paste("The number of component series in your data",
			  "is larger than their length. Are you sure",
			  "this is what you want?"))

 if (!is.numeric(y)) stop("Argument 'y' (data matrix) must be numeric.")

 if (n < 2) stop("Argument 'y' (data matrix) must contain at least two rows.")

 if (any(y == 0) && factors == 0) {
  myoffset <- sd(y)/10000
  warning(paste("Argument 'y' (data matrix) contains one or more zero(s) ",
	       	"and the number of factors is also zero. ",
		"I am adding an offset constant of size ", myoffset,
		" to do the auxiliary mixture sampling. If you want to ",
		"avoid this, you might consider de-meaning the returns ",
		"before calling this function.", sep=""))
 } else myoffset <- 0

 # Some error checking for draws
 if (!is.numeric(draws) | draws < 1) {
  stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
 } else {
  draws <- as.integer(draws)
 }

 # Some error checking for factors
 if (!is.numeric(factors) | factors < 0) {
  stop("Argument 'factors' (number of latent factors) must be a single number >= 0.")
 } else {
  factors <- as.integer(factors)
 }

 # Some error checking for burnin
 if (!is.numeric(burnin) | burnin < 0) {
  stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
 } else {
  burnin <- as.integer(burnin)
 }

 # Some error checking for zeromean/designmatrix
 if (isTRUE(zeromean)) {
   designmatrix <- NA
 } else if (identical(zeromean, FALSE)) {
   designmatrix <- "ar0"
 } else {
   stop("Argument 'zeromean' must be a single logical value.")
 }
 if (any(is.na(designmatrix))) {
   designmatrix <- NA
 } else if (identical(designmatrix, "ar0")) {
   ;
 } else {
   stop("Argument 'designmatrix' must be either NA or the character string 'ar0'.")
 }
 model_mean <- identical(designmatrix, "ar0")

 # Some error checking for priorh0idi
 if (length(priorh0idi) == 1) priorh0idi <- rep(priorh0idi, m)
 if (length(priorh0idi) != m) stop("Argument 'priorh0idi' must be of length 1 or ncol(y).")
 priorh0idi[remember <- priorh0idi == "stationary"] <- -1
 priorh0idi[!remember] <- as.numeric(priorh0idi[!remember])^2
 priorh0idi <- as.numeric(priorh0idi)
 if (any(priorh0idi[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")

 # Some error checking for priorh0fac
 if (length(priorh0fac) == 1) priorh0fac <- rep(priorh0fac, factors)
 if (length(priorh0fac) != factors) stop("Argument 'priorh0fac' must be of length 1 or factors.")
 priorh0fac[remember <- priorh0fac == "stationary"] <- -1
 priorh0fac[!remember] <- as.numeric(priorh0fac[!remember])^2
 priorh0fac <- as.numeric(priorh0fac)
 if (any(priorh0fac[!remember] < 0)) stop("Argument 'priorh0fac' must not contain negative values.")

 priorh0 <- c(priorh0idi, priorh0fac)


if (is.numeric(interweaving) && length(interweaving) == 1) {
 interweaving <- as.integer(interweaving)
} else {
 stop("Argument 'interweaving' must contain a single numeric value.")
}

if (interweaving != 0 & interweaving != 1 & interweaving != 2 & interweaving != 3 & interweaving != 4 & interweaving != 5 & interweaving != 6 & interweaving != 7) {
 stop("Argument 'interweaving' must be one of: 0, 1, 2, 3, 4.")
}


 # Some error checking for heteroskedastic
 if (length(heteroskedastic) == 1) heteroskedastic <- rep(heteroskedastic, m + factors)
 if (length(heteroskedastic) == 2) heteroskedastic <- c(rep(heteroskedastic[1], m), rep(heteroskedastic[2], factors))
 if (length(heteroskedastic) != m + factors) stop("Argument 'heteroskedastic' must be of length 1, 2, or (ncol(y) + factors).")
 if (!is.logical(heteroskedastic)) stop("Argument 'heteroskedastic' must be a vector containing only logical values.")
 if (is.null(heteroskedastic)) heteroskedastic <- rep(TRUE, m + factors)
 if (!all(heteroskedastic[m+seq_len(factors)])) {
  if (interweaving == 2L || interweaving == 4L) {
   warning("Cannot do deep interweaving if (some) factors are homoskedastic. Setting 'interweaving' to 3.")
   interweaving <- 3L
  }
 }

 if (!all(heteroskedastic)) {
  if (any(is.na(priorhomoskedastic))) {
   priorhomoskedastic <- matrix(c(1.1, 0.055), byrow = TRUE, nrow = m, ncol = 2)
   warning(paste0("Argument 'priorhomoskedastic' must be a matrix with dimension c(m, 2) if some of the
		  elements of 'heteroskedastic' are FALSE. Setting priorhomoskedastic to a matrix with
		  all rows equal to c(", priorhomoskedastic[1], ", ", priorhomoskedastic[2], ")."))
  }
  if (!is.matrix(priorhomoskedastic) || nrow(priorhomoskedastic) != m ||
      ncol(priorhomoskedastic) != 2 || any(priorhomoskedastic <= 0)) {
   stop("Argument 'priorhomoskedastic' must be a matrix with positive entries and dimension c(m, 2).")
  }
 }
 priorhomoskedastic <- as.matrix(priorhomoskedastic)

 # Some error checking for the prior parameters
 if (!is.numeric(priormu) | length(priormu) != 2) {
  stop("Argument 'priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
 }

 if (!is.numeric(priorphiidi) | length(priorphiidi) != 2) {
  stop("Argument 'priorphiidi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
 }

 if (!is.numeric(priorphifac) | length(priorphifac) != 2) {
  stop("Argument 'priorphifac' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
 }

 priorphi <- c(priorphiidi, priorphifac)

 if (!is.numeric(priorbeta) | length(priorbeta) != 2) {
  stop("Argument 'priorbeta' (mean and sd for the Gaussian prior for beta) must be numeric and of length 2.")
 }

 if (!is.numeric(priorsigmaidi) | any(priorsigmaidi <= 0)) {
  stop("Argument 'priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
 }

 if (length(priorsigmaidi) == 1) {
  priorsigmaidi <- rep(priorsigmaidi, m)
 } else if (length(priorsigmaidi) == m) {
  priorsigmaidi <- priorsigmaidi
 } else {
  stop("Argument 'priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must of length 1 or ncol(y).")
 }

 if (!is.numeric(priorsigmafac) | any(priorsigmafac <= 0)) {
  stop("Argument 'priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
 }

 if (length(priorsigmafac) == 1) {
  priorsigmafac <- rep(priorsigmafac, factors)
 } else if (length(priorsigmafac) == factors) {
  priorsigmafac <- priorsigmafac
 } else {
  stop("Argument 'priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must of length 1 or factors.")
 }

 priorsigma <- c(priorsigmaidi, priorsigmafac)

 if (!is.numeric(priorng) || length(priorng) != 2 || any(priorng <= 0)) {
  stop("Argument 'priorng' (prior hyperhyperparameters for Normal-Gamma prior) must be numeric and of length 2.")
 }

 cShrink <- priorng[1]
 dShrink <- priorng[2]

 if (!is.numeric(priorfacload) || any(priorfacload <= 0)) {
  stop("Argument 'priorfacload' must be numeric and positive.")
 }
 if (!(priorfacloadtype %in% c("normal", "rowwiseng", "colwiseng", "dl"))) {
  stop("Argument 'priorfacloadtype' must be one of: 'normal', 'rowwiseng',
       'colwiseng', 'dl'.")
 }

 if (!missing(startfacloadvar)) {
  if (!is.numeric(startfacloadvar) || !is.matrix(startfacloadvar) || nrow(startfacloadvar) != m ||
     ncol(startfacloadvar) != factors || any(startfacloadvar <= 0))
   stop("If argument 'startfacloadvar' is provided, it must be a matrix of appropriate dimensions with positive real entries.")
  if (priorfacloadtype == "normal")
    warning("Because priorfacloadtype is 'normal', the values passed via 'startfacloadvar' are being ignored.")
  startfacloadvarUsed <- startfacloadvar
 } else {
  startfacloadvarUsed <- matrix(1, nrow = m, ncol = factors)
 }

 if(is.matrix(priorfacload)) {
  if(nrow(priorfacload) != m || ncol(priorfacload) != factors) {
   stop("If argument 'priorfacload' is a matrix, it must be of appropriate dimensions.")
  }
  if (priorfacloadtype == "normal") {
   pfl <- 1L
   starttau2 <- priorfacload^2
   aShrink <- NA
   cShrink <- NA
   dShrink <- NA
  } else if (priorfacloadtype == "rowwiseng") {
   pfl <- 2L
   starttau2 <- startfacloadvarUsed
   aShrink <- priorfacload[,1]
   warning("Only first column of 'priorfacload' is used.'")
   cShrink <- rep(cShrink, m)
   dShrink <- rep(dShrink, m)
  } else if (priorfacloadtype == "colwiseng") {
   pfl <- 3L
   starttau2 <- startfacloadvarUsed
   aShrink <- priorfacload[1,]
   warning("Only first row of 'priorfacload' is used.'")
   cShrink <- rep(cShrink, factors)
   dShrink <- rep(dShrink, factors)
  } else if (priorfacloadtype == "dl") {
   pfl <- 4L
   starttau2 <- startfacloadvarUsed
   aShrink <- priorfacload[1,1]
   warning("Only first element of 'priorfacload' is used.'")
   cShrink <- NA
   dShrink <- NA
  }
 } else {
  if (length(priorfacload) != 1) {
   stop("If argument 'priorfacload' isn't a matrix, it must be a single number.")
  }
  if (priorfacloadtype == "normal") {
   pfl <- 1L
   starttau2 <- matrix(priorfacload^2, nrow = m, ncol = factors)
   aShrink <- NA
   cShrink <- NA
   dShrink <- NA
  } else if (priorfacloadtype == "rowwiseng") {
   pfl <- 2L
   starttau2 <- startfacloadvarUsed
   aShrink <- rep(priorfacload, m)
   cShrink <- rep(cShrink, m)
   dShrink <- rep(dShrink, m)
  } else if (priorfacloadtype == "colwiseng") {
   pfl <- 3L
   starttau2 <- startfacloadvarUsed
   aShrink <- rep(priorfacload, factors)
    cShrink <- rep(cShrink, factors)
    dShrink <- rep(dShrink, factors)
  } else if (priorfacloadtype == "dl") {
   pfl <- 4L
   starttau2 <- startfacloadvarUsed
   aShrink <- priorfacload
   cShrink <- NA
   dShrink <- NA
  }
 }

shrinkagepriors <- list(a = aShrink,
			c = cShrink,
			d = dShrink)

if(facloadtol < 0){
  stop("Argument 'facloadtol' (tolerance for the factor loadings) must be >=0.")
}

 # Some error checking for thin
 if (!is.numeric(thin) | thin < 1) {
  stop("Argument 'thin' (thinning parameter for the latent log variances and parameters) must be a single number >= 1.")
 } else {
  thin <- as.integer(thin)
 }

 # Some error checking for keeptime

 if (is.numeric(keeptime)) {
  thintime <- as.integer(keeptime)
 } else {
  thintime <- switch(keeptime,
                all = 1L,
                last = -1L)
 }
 if (is.null(thintime)) stop("Argument 'keeptime' must be either numeric > 0, \"all\" or \"last\".")

 # Some error checking for expert
 if (missing(expert)) {
  para <- 3L ; parameterization <- 'GIS_C'
  mhcontrol <- -1
  gammaprior <- TRUE
  truncnormal <- FALSE
  mhsteps <- 2L
  B011 <- 10^8
  B022 <- 10^12
 } else {
  expertnames <- names(expert)
  if (!is.list(expert) | is.null(expertnames) | any(expertnames == ""))
   stop("Argument 'expert' must be a named list with nonempty names.")
  if (length(unique(expertnames)) != length(expertnames))
   stop("No duplicate elements allowed in argument 'expert'.")
  allowednames <- c("parameterization", "mhcontrol", "gammaprior",
		    "truncnormal", "mhsteps", "proposalvar4sigmaphi",
		    "proposalvar4sigmatheta")
  exist <- pmatch(expertnames, allowednames)
  if (any(is.na(exist)))
   stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))

  expertenv <- list2env(expert)

  if (exists("parameterization", expertenv)) {
   parameterization <- expert[["parameterization"]]
   if (!is.character(parameterization) | is.na(parameterization)) {
    stop("Argument 'parameterization' must be either 'centered', 'noncentered', 'GIS_C', or 'GIS_NC'.")
   }
   switch(parameterization,
    centered = para <- 1L,
    noncentered = para <- 2L,
    GIS_C = para <- 3L,
    GIS_NC = para <- 4L,
    stop("Unknown parameterization. Currently you can only use 'centered', 'noncentered', 'GIS_C', and 'GIS_NC'.")
   )
  } else {
   para <- 3L ; parameterization <- 'GIS_C'
  }

  # Remark: mhcontrol < 0 means independence proposal,
  #         mhcontrol > 0 controls stepsize of log-random-walk proposal
  if (exists("mhcontrol", expertenv)) {
   mhcontrol <- expert[["mhcontrol"]]
   if (!is.numeric(mhcontrol) | length(mhcontrol) != 1)
    stop("Argument 'mhcontrol' must be a single number.")
  } else {
   mhcontrol <- -1
  }

  # use a Gamma prior for sigma^2 in C?
  if (exists("gammaprior", expertenv)) {
   gammaprior <- expert[["gammaprior"]]
   if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
  } else {
   gammaprior <- TRUE
  }

  # use a truncated normal as proposal? (or normal with rejection step)
  if (exists("truncnormal", expertenv)) {
   truncnormal <- expert[["truncnormal"]]
   if (!is.logical(truncnormal)) stop("Argument 'truncnormal' must be TRUE or FALSE.")
  } else {
   truncnormal <- FALSE
  }

  if (exists("mhsteps", expertenv)) {
   mhsteps <- as.integer(expert[["mhsteps"]])
   if (mhsteps != 2L & mhsteps != 1L & mhsteps != 3L) stop("mhsteps must be 1, 2, or 3")
   if (mhsteps != 2L & mhcontrol >= 0)
    stop("Log normal random walk proposal currently only implemented for mhsteps==2.")
   if (mhsteps != 2L & !isTRUE(gammaprior))
    stop("Inverse Gamma prior currently only implemented for mhsteps==2.")
  } else {
   mhsteps <- 2L
  }

  # prior for ridge _proposal_ (variance of sigma*phi)
  if (exists("proposalvar4sigmaphi", expertenv)) {
   B011 <- expert[["proposalvar4sigmaphi"]]
   if (!is.numeric(B011) | length(B011) != 1 | B011 <= 0)
    stop("Argument 'proposalvar4sigmaphi' must be a positive number.")
  } else {
   B011 <- 10^8
  }

  # prior for ridge _proposal_ (variance of sigma*theta)
  if (exists("proposalvar4sigmatheta", expertenv)) {
   B022 <- expert[["proposalvar4sigmatheta"]]
   if (!is.numeric(B022) | length(B022) != 1 | B022 <= 0)
    stop("Argument 'proposalvar4sigmatheta' must be a positive number.")
  } else {
   B022 <- 10^12
  }
 }

# Some input checking for startfacload
if (missing(startfacload)) {
 startfacload <- matrix(rnorm(m*factors, sd = .5)^2, nrow=m, ncol=factors)
# if (factors >= 1) for (i in 1:factors) startfacload[i,] <- c(rep(1/i,i) + rnorm(i, sd=.1), rep(0,factors-i))
# startfacload[-(1:factors),] <- 1/factors + rnorm((m-factors)*factors, sd=.1)
} else {
 if (factors == 1L && is.numeric(startfacload) && is.vector(startfacload) && length(startfacload) == m)
  startfacload <- matrix(startfacload, ncol = 1L)
 if (!is.numeric(startfacload) || !is.matrix(startfacload) ||
     (nrow(startfacload) != m && factors >= 1) || ncol(startfacload) != factors)
  stop("Argument 'startfacload' must be a numeric matrix of dimension c(ncol(y), factors).")
}

 # Some input checking for startpara
 if (missing(startpara)) {
  startpara <- list(mu = c(rep(-3, m) + rnorm(m), rep(0, factors)),
		    phi = c(rep(.8, m), rep(.8, factors)) + pmin(rnorm(m + factors, sd=.06), .095),
		    sigma = rep(.1, m + factors) + rgamma(m + factors, 1, 10))
 } else {
  if (!is.numeric(startpara) || !is.matrix(startpara) || nrow(startpara) != 3 || ncol(startpara) != m + factors)
   stop("Argument 'startpara' must be a numeric matrix of dimension c(3, ncol(y) + factors).")

  if (any(startpara[1, m + seq_len(factors)] != 0)) {
   warning("Some of the levels of the factor log variance starting values are different from zero. Setting those to zero for you.")
   startpara[1, m + seq_len(factors)] <- 0
  }

  if (any(abs(startpara[2,]) >= 1))
   stop("All elements of 'startpara[2,]' must be between -1 and 1.")

  if (any(startpara[3,] <= 0))
   stop("All elements of 'startpara[3,]' must be positive.")

  startpara <- list(mu = startpara[1,], phi = startpara[2,], sigma = startpara[3,]) # this is what the sampler expects
 }

 # Some input checking for startlogvar
 if (missing(startlogvar)) {
  startlogvar <- matrix(startpara[["mu"]][1] + rnorm(n*(m + factors)), n, m + factors)
 } else {
  if (!is.numeric(startlogvar) || !is.matrix(startlogvar) ||
      nrow(startlogvar) != nrow(y) || ncol(startlogvar) != factors + m)
   stop("Argument 'startlogvar' must be a numeric matrix of dimension c(nrow(y), ncol(y) + factors).")
 }

 if (missing(startlogvar0)) {
  startlogvar0 <- startpara[["mu"]][1] + rnorm(m + factors)
 } else {
  if (!is.numeric(startlogvar0) || length(startlogvar0) != (m + factors))
   stop("Argument 'startlogvar0' must be a vector of length ncol(y) + factors.")
 }

# Some input checking for startfac
if (missing(startfac)) {
 startfac <- matrix(rnorm(factors*n, 0, sd=.1), nrow=factors)
} else {
 if (factors == 1L && is.numeric(startfac) && is.vector(startfac) && length(startfac) == n)
  startfac <- matrix(startfac, nrow = 1L)
 if (!is.numeric(startfac) || !is.matrix(startfac) ||
     nrow(startfac) != factors || (ncol(startfac) != n && factors >= 1))
  stop("Argument 'startfac' must be a numeric matrix of dimension c(factors, nrow(y)).")
}

if (is.numeric(runningstore) && length(runningstore) == 1 && runningstore >= 0 && runningstore <= 6) {
 runningstore <- as.integer(runningstore)
} else {
 stop("Argument 'runningstore' must be one of: 0, 1, 2, 3, 4, 5, 6.")
}

if (runningstore > 1 && factors == 0) {
 runningstore <- 1L
 warning("Setting 'runningstore' to something greater than 1 isn't supported in a no-factor model. I set 'runningstore' to 1 for you.")
}

if (!is.numeric(runningstorethin) || length(runningstorethin) != 1 || runningstorethin < 1) {
 stop("Argument 'runningstorethin' must be a positive integer.")
 runningstorethin <- as.integer(runningstorethin)
}

if (!is.numeric(runningstoremoments) || length(runningstoremoments) != 1 || runningstoremoments < 1 || runningstoremoments > 4) {
 stop("Argument 'runningstoremoments' must be one of: 1, 2, 3, 4.")
 runningstoremoments <- as.integer(runningstoremoments)
}


 if (!quiet) {
  cat(paste("\nCalling ", factors, "-factor MCMC sampler ",
#	    "with ", draws+burnin, " iter ",
	    "for ", ncol(y), " series of length ",
	    nrow(y), ".\n\n", sep = ""))
  flush.console()
 }

 ## Hack to prevent console flushing problems with Windows
 #if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet
 myquiet <- quiet


if (is.matrix(restrict)) {
 if (any(dim(restrict) != c(m, factors)) || any(is.na(restrict)) || !is.logical(restrict))
     stop("Argument 'restrict' must be an appropriate logical matrix or \"none\"/\"upper\"/\"auto\".")
 if (any(rowSums(restrict) == factors))
     stop("Argument 'restrict' can't have rows where all elements are TRUE.")
 restr <- matrix(FALSE, nrow = m, ncol = factors)
 restr[restrict == TRUE] <- TRUE
}

if (is.character(restrict)) {
 if (length(restrict) > 1) stop("Argument 'restrict' must be an appropriate matrix or \"none\"/\"upper\"/\"auto\".")
 restr <- matrix(FALSE, nrow = m, ncol = factors)
 if (restrict == "upper") restr[upper.tri(restr)] <- TRUE
 if (restrict == "auto") restr <- findrestrict(y, factors = factors)
}

if (interweaving %in% c(1, 2) && any(diag(restr) == TRUE)) {
  stop("Setting 'interweaving' to either 1 or 2 and restricting the diagonal elements of the factor loading matrix are not allowed at the same time.")
}

if (length(samplefac) != 1L || !is.logical(samplefac)) {
  stop("Argument 'samplefac' must be logical and of length one.")
}

if (!isTRUE(samplefac) && isTRUE(signswitch)) {
  signswitch <- FALSE
  warning("Turning sign switching off because 'samplefac' is set to FALSE.")
}

if (!isTRUE(samplefac) && interweaving != 0L) {
  interweaving <- 0L
  warning("Turning interweaving off because 'samplefac' is set to FALSE.")
}

if (!isTRUE(samplefac) && isTRUE(signident)) {
  signident <- FALSE
  warning("Turning ex-post sign-identification off because 'samplefac' is set to FALSE.")
}

restrinv <- matrix(as.integer(!restr), nrow = nrow(restr), ncol = ncol(restr))

startval <- list(facload = startfacload,
		 fac = startfac,
		 para = startpara,
		 latent = startlogvar,
		 latent0 = startlogvar0,
		 tau2 = starttau2)

auxstore <- FALSE

res <- .Call("sampler", t(y), draws, burnin, startval,
             priormu[1], priormu[2]^2, priorphi, priorsigma,
             priorbeta, model_mean, shrinkagepriors,
             thin, auxstore, thintime, myquiet, para,
             mhsteps, B011, B022, mhcontrol, gammaprior,
             myoffset, truncnormal,
             restrinv, interweaving, signswitch, runningstore,
             runningstorethin, runningstoremoments, pfl,
             heteroskedastic, priorhomoskedastic, priorh0, samplefac, facloadtol,
             PACKAGE = "factorstochvol")

res$y <- y

res$config <- list(draws = draws, burnin = burnin, thin = thin,
		    keeptime = keeptime, runningstore = runningstore,
		    runningstorethin = runningstorethin,
		    runningstoremoments = runningstoremoments,
		    restrict = restr, interweaving = interweaving,
		    signswitch = signswitch,
		    heteroskedastic = heteroskedastic,
        designmatrix = designmatrix,
		    expert = list(parameterization = parameterization,
				  mhcontrol = mhcontrol,
				  gammaprior = gammaprior,
				  truncnormal = truncnormal,
				  mhsteps = mhsteps,
				  B011 = B011, B022 = B022),
		    startpara = startpara,
		    startlogvar = startlogvar,
		    startlogvar0 = startlogvar0,
		    startfacload = startfacload,
		    startfac = startfac,
		    startfacloadvar = starttau2,
        samplefac = samplefac)
 res$priors <- list(priormu = priormu,
		    priorphiidi = priorphiidi,
		    priorphifac = priorphifac,
		    priorsigmaidi = priorsigmaidi,
		    priorsigmafac = priorsigmafac,
		    priorfacload = priorfacload,
		    priorfacloadtype = priorfacloadtype,
		    priorng = priorng, priorbeta = priorbeta,
		    priorh0idi = priorh0idi, priorh0fac = priorh0fac,
		    priorhomoskedastic = priorhomoskedastic)

 if (!auxstore) {
  res$tau2 <- NULL
  res$lambda2 <- NULL
  res$mixind <- NULL
 }

 if (NROW(res$beta) == 0) {
   res$beta <- NULL
 } else {
   rownames(res$beta) <- paste0("mean_", seq_len(NROW(res$beta)))
 }

 dimnames(res$para) <- list(c("mu", "phi", "sigma"), NULL, NULL)

 if (!quiet) {
  cat("\n\nReorganizing runningstores... ")
  flush.console()
 }

 mynames <- c("mean", "m2", "m3", "m4")[seq_len(runningstoremoments)]

 if (runningstore >= 6L) {  # stored running communalities, correlations, covariances, factor sds, factors, and factor log variances
   res$runningstore$com <- standardizer(force(array(unlist(res$runningstore$com, FALSE, FALSE),
     dim = c(dim(res$runningstore$com$mean), runningstoremoments),
     dimnames = list(NULL, NULL, mynames))))
 } else {
   res$runningstore$com <- NULL
 }

 if (runningstore >= 5L) {  # stored running correlations, covariances, factor sds, factors, and factor log variances
   tmpnames <- apply(expand.grid(1:m, 1:m)[lower.tri(diag(m), diag = FALSE),], 1, paste, collapse = "_")
   res$runningstore$cor <- standardizer(force(array(unlist(res$runningstore$cor, FALSE, FALSE),
     dim = c(dim(res$runningstore$cor$mean), runningstoremoments),
     dimnames = list(NULL, tmpnames, mynames))))
 } else {
   res$runningstore$cor <- NULL
 }

 if (runningstore >= 4L) {  # stored running covariances (and sqrt(diag(covariances)),
                            # factor sds, factors, and factor log variances
   tmpnames <- apply(expand.grid(1:m, 1:m)[lower.tri(diag(m), diag = TRUE),], 1, paste, collapse = "_")

   res$runningstore$cov <- standardizer(force(array(unlist(res$runningstore$cov, FALSE, FALSE),
     dim = c(dim(res$runningstore$cov$mean), runningstoremoments),
     dimnames = list(NULL, tmpnames, mynames))))

   res$runningstore$vol <- standardizer(force(array(unlist(res$runningstore$vol, FALSE, FALSE),
     dim = c(dim(res$runningstore$vol$mean), runningstoremoments),
     dimnames = list(NULL, NULL, mynames))))
 } else {
   res$runningstore$cov <- NULL
   res$runningstore$vol <- NULL
 }

 if (runningstore >= 3L) {  # stored running factor sds, factors, and factor log variances
   res$runningstore$sd <- standardizer(force(array(unlist(res$runningstore$sd, FALSE, FALSE),
     dim = c(dim(res$runningstore$logvar$mean), runningstoremoments),
     dimnames = list(NULL, NULL, mynames))))
 } else {
   res$runningstore$sd <- NULL
 }

 if (runningstore >= 2L) {  # stored running factors and factor log variances
   res$runningstore$fac <- standardizer(force(array(unlist(res$runningstore$fac, FALSE, FALSE),
     dim = c(dim(res$runningstore$fac$mean), runningstoremoments),
     dimnames = list(NULL, NULL, mynames))))
 } else {
   res$runningstore$fac <- NULL
 }

 if (runningstore >= 1L) {  # stored running factor log variances
    res$runningstore$logvar <- standardizer(force(array(unlist(res$runningstore$logvar, FALSE, FALSE),
      dim = c(dim(res$runningstore$logvar$mean), runningstoremoments),
      dimnames = list(NULL, NULL, mynames))))
 } else {
  res$runningstore <- NULL
 }

 if (!quiet) {
  cat("Done!\n")
  flush.console()
 }

 class(res) <- "fsvdraws"

 if (signident) {
 if (!quiet) {
  cat("Ex-post sign-identification... ")
  flush.console()
 }

 res <- signident(res)
 }

 if (!quiet) {
  cat("Done!\n")
  flush.console()
 }
 res
}


#' Predicts means and variances conditionally on the factors
#'
#' \code{predcond} simulates from the posterior predictive distribution
#' of the data, conditionally on realized values of the factors. This
#' has the advantage that the predictive density can be written as
#' the product of the marginals but introduces sampling uncertainty
#' that grows with the number of factors used.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param ahead Vector of timepoints, indicating how many steps
#' to predict ahead.
#' @param each Single integer (or coercible to such) indicating how
#' often should be drawn from the posterior predictive distribution
#' for each draw that has been stored during MCMC sampling.
#' @param ... Ignored.
#'
#' @return List of class \code{fsvpredcond} containing two elements:
#' \itemize{
#' \item{means}{Array containing the draws of the predictive means.}
#' \item{vols}{Array containing the draws of the predictive volatilities (square root of variances).}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- fsvsim(n = 500, series = 4, factors = 1) # simulate
#' res <- fsvsample(sim$y, factors = 1) # estimate
#'
#' # Predict 1 day ahead:
#' predobj <- predcond(res, each = 5)
#'
#' # Draw from the predictive distribution:
#' preddraws <- matrix(rnorm(length(predobj$means[,,1]),
#'                     mean = predobj$means[,,1],
#'                     sd = predobj$vols[,,1]), nrow = 4)
#'
#' # Visualize the predictive distribution
#' pairs(t(preddraws), col = rgb(0,0,0,.1), pch = 16)
#' }
#'
#' @family predictors
#'
#' @export
predcond <- function(x, ahead = 1, each = 1, ...) {

 if (!is(x, "fsvdraws")) stop("Argument 'x' must be of class 'fsvdraws'.")

 if (!is.vector(ahead) | !is.numeric(ahead) | any(is.na(ahead)))
  stop("Argument 'ahead' must be a numeric vector, NAs are not allowed.")

 ahead <- as.integer(ahead)
 ahead <- sort(ahead)

 if (any(ahead < 1)) stop("All elements of 'ahead' must be greater or equal to 1.")

 if (!is.vector(each) || !is.numeric(each) || any(is.na(each)) || length(each) != 1 || each < 1)
  stop("Argument 'each' must be a single number >= 1, NAs are not allowed.")

 each <- as.integer(each)

 res <- .Call("predict", x, ahead, each, PACKAGE = "factorstochvol")

 if (!is.null(x$beta)) {
  for (j in seq_len(NROW(x$beta))) {
   # res$beta[j, ] is one distribution so it doesn't even matter how it is recycled
   res$means[j, , ] <- res$means[j, , ] + x$beta[j, ]
  }
 }

 dimnames(res$means) <- list(NULL, NULL, ahead = ahead)
 dimnames(res$vols) <- list(NULL, NULL, ahead = ahead)

 class(res) <- c("fsvpredcond")
 res
}


myrgig <- function(n = 1, lambda, chi, psi) {

  ## --------------------------------------------------------------------
  ## Generate GIG distributed variables (imported from package GIGrvg)
  ##
  ## density proportional to
  ##    f(x) = x^{lambda-1} e^{-1/2 (chi/x+psi x)}
  ##
  ##       x >= 0
  ## --------------------------------------------------------------------
  ## Arguments:
  ##
  ##   n ....... sample size
  ##   lambda .. parameter for distribution
  ##   chi   ... parameter for distribution
  ##   psi   ... parameter for distribution
  ## --------------------------------------------------------------------

  ## generate sample
  .Call("my_rgig", n, lambda, chi, psi, PACKAGE = "factorstochvol")
}

standardizer <- function(obj) {
 if (dim(obj)[3] >= 2L) {
  tmpmean <- obj[,,"mean",drop=FALSE]
  tmpsd <- sqrt(obj[,,"m2",drop=FALSE] - tmpmean^2)
  obj[,,"m2"] <- tmpsd
  dimnames(obj)[[3]][dimnames(obj)[[3]] == "m2"] <- "sd"
 }
 if (dim(obj)[3] >= 3L) {
  tmpskew <- (obj[,,"m3",drop=FALSE] - 3*tmpmean*tmpsd^2 - tmpmean^3) / tmpsd^3
  obj[,,"m3"] <- tmpskew
  dimnames(obj)[[3]][dimnames(obj)[[3]] == "m3"] <- "skew"
 }
 if (dim(obj)[3] >= 4L) {
  tmpkurt <- (obj[,,"m4",drop=FALSE] - tmpmean * (4*tmpsd^3*tmpskew + 6*tmpmean*tmpsd^2 + tmpmean^3)) / tmpsd^4
  obj[,,"m4"] <- tmpkurt
  dimnames(obj)[[3]][dimnames(obj)[[3]] == "m4"] <- "kurt"
 }
 # sanity check:
 #stopifnot(all.equal(dimnames(obj)[[3]], c("min", "mean", "m2", "m3", "m4", "max")))
 #dimnames(obj)[[3]] <- c("min", "mean", "sd", "skew", "exkurt", "max")
 obj
}

## R wrapper for dmvnorm ("vectorized")
## Caveat! Does not do proper input checking (i.e. positive definiteness, etc.)
## No recycling!

vecdmvnorm <- function(x, means, vars, log = FALSE, tol = 10^6 * .Machine$double.eps) {
 if (!is.matrix(x) || !is.numeric(x)) stop("Argument 'x' must be a numeric matrix.")
 if (!is.matrix(means) || !is.numeric(means)) stop("Argument 'means' must be a numeric matrix.")
 if (!is.array(vars) || !is.numeric(vars)) stop("Argument 'vars' must be a numeric array.")
 if (!isSymmetric(vars[,,1], tol)) stop("Each slice of 'vars' must be symmetric.")
 if (!all.equal(dim(x), dim(means))) stop("Dimensions of 'x' and 'means' do not match.")
 m <- nrow(x)
 n <- ncol(x)
 if (nrow(vars) != m || ncol(vars) != m) stop("Number of rows and columns of 'vars' must be equal to nrow(x).")
 if (dim(vars)[3] != n) stop("Argument 'vars' must have exactly ncol(x) slices.")

 .Call("dmvnorm", x, means, vars, log, PACKAGE = "factorstochvol")
}
