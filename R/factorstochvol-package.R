#' Bayesian Estimation of (Sparse) Latent Factor Stochastic
#' Volatility Models through MCMC
#'
#'This packages provides a Markov chain Monte Carlo (MCMC) sampler
#' for fully Bayesian estimation of latent factor stochastic volatility
#' models. Sparsity can be achieved through the usage of Normal-Gamma
#' priors on the factor loadings matrix.
#'
#' In recent years, multivariate factor stochastic volatility (SV)
#' models have been increasingly used to analyze financial and economic
#' time series because they can capture joint (co-)volatility dynamics
#' by a small number of latent time-varying factors. The main advantage
#' of such a model is its parsimony, as all variances and covariances
#' of a time series vector are governed by a low-dimensional common factor
#' with the components following independent SV models. For problems of
#' this kind, MCMC is a very efficient estimation method, it is however
#' associated with a considerable computational burden when the number
#' of assets is moderate to large. To overcome this, the latent volatility
#' states are drawn "all without a loop" (AWOL), ancillarity-sufficiency
#' interweaving strategies (ASIS) are applied to sample the univariate
#' components as well as the factor loadings. Thus, this package can
#' be applied directly estimate time-varying covariance and correlation
#' matrices for medium-and high-dimensional time series. To guarantee
#' sparsity, a hierarchical Normal-Gamma prior can be used for the
#' factor loadings matrix which shrinks the unnecessary factor loadings
#' towards zero. 
#'
#' @note This package is currently in active development; the interface
#' of some of the functions might change.
#' Moreover, even though I tried to carefully check everything,
#' factorstochvol may still contain
#' typos, inconsistencies, or even bugs. Your comments and suggestions
#' are warmly welcome!
#'
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @references
#' Kastner, G., Frühwirth-Schnatter, S., Lopes, H. F. (2016). Efficient
#' Bayesian inference for multivariate factor stochastic volatility models.
#' \emph{Report 128, Research Report Series of the Institute of Statistics
#' and Mathematics, WU Vienna University of Economics and Business},
#' \url{http://epub.wu.ac.at/4875/}
#' 
#' Kastner, G. (2016). Sparse Bayesian time-varying covariance estimation
#' in many dimensions. \emph{Working Paper}
#' 
#' Kastner, G. and Frühwirth-Schnatter, S. (2014). Ancillarity-sufficiency
#' interweaving strategy (ASIS) for boosting MCMC estimation of stochastic
#' volatility models. \emph{Computational Statistics and Data Analysis},
#' \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.
#'
#' @keywords package models ts
#'
#' @seealso \code{\link[stochvol:stochvol-package]{stochvol}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' 
#' # simulate data from a (small) factor SV model:
#' sim <- fsvsim(series = 5, factors = 2)
#' 
#' # estimate the model (CAVEAT: only few draws!)
#' res <- fsvsample(sim$y, factors = 2, draws = 2000,
#'                  burnin = 500, runningstore = 6)
#' 
#' # plot implied volas overtime:
#' voltimeplot(res)
#' 
#' # plot correlation matrix at some points in time:
#' par(mfrow = c(2,2))
#' corimageplot(res, seq(1, nrow(sim$y), length.out = 4),
#'              fsvsimobj = sim, plotCI = 'circle',
#'              plotdatedist = -2)
#' 
#' 
#' # plot (certain) covariances and correlations over time
#' par(mfrow = c(2,1))
#' covtimeplot(res, 1)
#' cortimeplot(res, 1)
#' 
#' # plot (all) correlations over time
#' corplot(res, fsvsimobj = sim, these = 1:10)
#' 
#' # plot factor loadings
#' par(mfrow = c(1,1))
#' facloadpointplot(res, fsvsimobj = sim)
#' facloadpairplot(res)
#' facloadcredplot(res)
#' facloaddensplot(res, fsvsimobj = sim)
#' 
#' # plot latent log variances
#' logvartimeplot(res, fsvsimobj = sim, show = "fac")
#' logvartimeplot(res, fsvsimobj = sim, show = "idi")
#' 
#' # plot communalities over time
#' comtimeplot(res, fsvsimobj = sim, show = 'joint')
#' comtimeplot(res, fsvsimobj = sim, show = 'series')
#' }
#'
#' @docType package
#' @name factorstochvol-package

NULL
