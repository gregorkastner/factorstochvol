/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2016-2021
 *     Darjus Hosszejni Copyright (C) 2019-2021
 *     Luis Gruber Copyright (C) 2021
 *
 *  This file is part of the R package factorstochvol: Bayesian Estimation
 *  of (Sparse) Latent Factor Stochastic Volatility Models
 *
 *  The R package factorstochvol is free software: you can redistribute
 *  it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or any
 *  later version of the License.
 *
 *  The R package factorstochvol is distributed in the hope that it will
 *  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the R package factorstochvol. If that is not the case,
 *  please refer to <http://www.gnu.org/licenses/>.
 */

#ifndef _UPDATE_FSV_H
#define _UPDATE_FSV_H

#include <RcppArmadillo.h>
#include <stochvol.h>

void update_fsv(arma::mat& facload,
                arma::mat& fac,
                arma::mat& logvar,
                arma::vec& logvar0,
                arma::mat& svpara,
                arma::mat& tau2,
                arma::vec& lambda2, 
                arma::umat& curmixind, 
                const arma::mat& y,
                const double& facloadtol,
                const arma::imat& restriction,
                const arma::uvec& facloadtunrestrictedelements,
                const arma::irowvec& nonzerospercol,
                const arma::icolvec& nonzerosperrow,
                const Rcpp::NumericVector priorh0,
                const bool& ngprior,
                const bool& columnwise,
                const Rcpp::NumericVector aShrink,
                const Rcpp::NumericVector cShrink,
                const Rcpp::NumericVector dShrink,
                const Rcpp::NumericMatrix priorhomoskedastic,
                const double& offset,
                const Rcpp::NumericVector heteroskedastic, 
                const int& interweaving,
                const stochvol::ExpertSpec_FastSV& expert_idi,
                const stochvol::ExpertSpec_FastSV& expert_fac,
                const std::vector<stochvol::PriorSpec>& prior_specs,
                const double& B011inv,
                const bool& samplefac,
                const bool& signswitch,
                const int& i
);

#endif
