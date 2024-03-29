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

#' @importFrom grDevices colorRampPalette rainbow rgb col2rgb palette
#' @importFrom graphics abline axis barplot image layout legend lines mtext par plot plot.new points symbols text title
#' @importFrom stats cov2cor density factanal median quantile rgamma rnorm runif rt sd ts.plot
#' @importFrom utils combn flush.console
#' @importFrom corrplot corrplot corrMatOrder
#' @importFrom Rcpp evalCpp
#' @importFrom stochvol logret paratraceplot svsample validate_and_process_expert
#' @importFrom GIGrvg rgig
NULL

