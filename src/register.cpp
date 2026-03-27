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

#include "update_fsv.h"
#include "../inst/include/factorstochvol.h"

using namespace Rcpp;

// validate (ensure exported C++ functions exist before calling them)
static int _factorstochvol_RcppExport_validate(const char* sig) { 
  static std::set<std::string> signatures;
  if (signatures.empty()) {
    signatures.insert("void(*update_fsv)(arma::mat&, arma::mat&, arma::mat&, arma::vec&, arma::mat&, arma::mat&, arma::vec&, arma::umat&, const arma::mat&, const double&, const arma::imat&, const arma::uvec&, const arma::irowvec&, const arma::icolvec&, const Rcpp::NumericVector, const bool&, const bool&, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericVector, const Rcpp::NumericMatrix, const double&, const Rcpp::NumericVector, const int&, const stochvol::ExpertSpec_FastSV&, const stochvol::ExpertSpec_FastSV&, const std::vector<stochvol::PriorSpec>&, const double&, const bool&, const bool&, const int&)");
  }
  return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _factorstochvol_RcppExport_registerCCallable() { 
  R_RegisterCCallable("factorstochvol", "update_fsv", (DL_FUNC)update_fsv);
  R_RegisterCCallable("factorstochvol", "_factorstochvol_RcppExport_validate", (DL_FUNC)_factorstochvol_RcppExport_validate);
  return R_NilValue;
}

