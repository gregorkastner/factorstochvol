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

#' Bayesian Estimation of (Sparse) Latent Factor Stochastic
#' Volatility Models through MCMC
#'
#' This packages provides a Markov chain Monte Carlo (MCMC) sampler
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
#'
#' @references
#' Kastner, G., Frühwirth-Schnatter, S., and Lopes, H.F. (2017).
#' Efficient Bayesian Inference for Multivariate Factor Stochastic Volatility Models.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26}(4), 905--917,
#' \doi{10.1080/10618600.2017.1322091}.
#'
#' Kastner, G. (2019). Sparse Bayesian Time-Varying Covariance Estimation
#' in Many Dimensions. \emph{Journal of Econometrics}, \bold{210}(1), 98--115.
#' \doi{10.1016/j.jeconom.2018.11.007}.
#'
#' Kastner, G. and Frühwirth-Schnatter, S. (2014). Ancillarity-Sufficiency
#' Interweaving Strategy (ASIS) for Boosting MCMC Estimation of Stochastic
#' Volatility Models. \emph{Computational Statistics and Data Analysis},
#' \doi{10.1016/j.csda.2013.01.002}.
#'
#' @keywords package models ts
#'
#' @seealso \code{\link[stochvol:stochvol-package]{stochvol}}
#'
#' @useDynLib factorstochvol, .registration = TRUE
#'
#' @examples
#' \donttest{
#' set.seed(1)
#'
#' # simulate data from a (small) factor SV model:
#' sim <- fsvsim(series = 5, factors = 2)
#'
#' # estimate the model (CAVEAT: only few draws!)
#' res <- fsvsample(sim$y, factors = 2, draws = 2000, burnin = 500)
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
