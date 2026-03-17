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

#ifndef _FACTORSTOCHVOL_H
#define _FACTORSTOCHVOL_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stochvol.h>

namespace factorstochvol {

using namespace Rcpp;

namespace {
void validateSignature(const char* sig) {
  Rcpp::Function require = Rcpp::Environment::base_env()["require"];
  require("factorstochvol", Rcpp::Named("quietly") = true);
  typedef int(*Ptr_validate)(const char*);
  static Ptr_validate p_validate = (Ptr_validate)
    R_GetCCallable("factorstochvol", "_factorstochvol_RcppExport_validate");
  if (!p_validate(sig)) {
    throw Rcpp::function_not_exported(
        "C++ function with signature '" + std::string(sig) + "' not found in factorstochvol");
  }
}
}

inline 
void update_fsv(arma::mat& facload, arma::mat& fac, arma::mat& logvar, arma::vec& logvar0, arma::mat& svpara, arma::mat& tau2, arma::vec& lambda2, arma::umat& curmixind, const arma::mat& y, const double& facloadtol, const arma::imat& restriction, const arma::uvec& facloadtunrestrictedelements, const arma::irowvec& nonzerospercol, const arma::icolvec& nonzerosperrow, const Rcpp::NumericVector priorh0, const bool& ngprior, const bool& columnwise, const Rcpp::NumericVector aShrink, const Rcpp::NumericVector cShrink, const Rcpp::NumericVector dShrink, const Rcpp::NumericMatrix priorhomoskedastic, const double& offset, const Rcpp::NumericVector heteroskedastic, const int& interweaving, const stochvol::ExpertSpec_FastSV& expert_idi, const stochvol::ExpertSpec_FastSV& expert_fac, const std::vector<stochvol::PriorSpec>& prior_specs, const double& B011inv, const bool& samplefac, const bool& signswitch, const int& i) {
  typedef void(*Ptr_update_fsv)(arma::mat&, arma::mat&, arma::mat&, arma::vec&, arma::mat&, arma::mat&, arma::vec&, arma::umat&, const arma::mat&, const double&, const arma::imat&, const arma::uvec&, const arma::irowvec&, const arma::icolvec&, const Rcpp::NumericVector, const bool&, const bool&, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericMatrix, const double&, const Rcpp::NumericVector, const int&, const stochvol::ExpertSpec_FastSV&, const stochvol::ExpertSpec_FastSV&, const std::vector<stochvol::PriorSpec>&, const double&, const bool&, const bool&, const int&);
  static Ptr_update_fsv p_update_fsv = NULL;
  if (p_update_fsv == NULL) {
    validateSignature("void(*update_fsv)(arma::mat&, arma::mat&, arma::mat&, arma::vec&, arma::mat&, arma::mat&, arma::vec&, arma::umat&, const arma::mat&, const double&, const arma::imat&, const arma::uvec&, const arma::irowvec&, const arma::icolvec&, const Rcpp::NumericVector, const bool&, const bool&, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericMatrix, const double&, const Rcpp::NumericVector, const int&, const stochvol::ExpertSpec_FastSV&, const stochvol::ExpertSpec_FastSV&, const std::vector<stochvol::PriorSpec>&, const double&, const bool&, const bool&, const int&)"); 
    p_update_fsv = (Ptr_update_fsv)R_GetCCallable("factorstochvol", "update_fsv");
  }
  {
    p_update_fsv(facload, fac, logvar, logvar0, svpara, tau2, lambda2, curmixind, y, facloadtol, restriction, facloadtunrestrictedelements, nonzerospercol, nonzerosperrow, priorh0, ngprior, columnwise, aShrink, cShrink, dShrink, priorhomoskedastic, offset, heteroskedastic, interweaving, expert_idi, expert_fac, prior_specs, B011inv, samplefac, signswitch, i);
  }
}

}

#endif