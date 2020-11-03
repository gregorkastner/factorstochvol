/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2016-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
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

#ifndef _PROGUTILS_H
#define _PROGUTILS_H

#include <RcppArmadillo.h>

double logdnormquot(double x, double y, double mu, double sigma);
double logspecialquot(double x, double y, double alpha, double beta, double c);

void store(const Rcpp::NumericMatrix &curfacload, Rcpp::NumericVector &facload,
           const Rcpp::NumericMatrix &curf,       Rcpp::NumericVector &f,
           const Rcpp::NumericMatrix &curh,       Rcpp::NumericVector &h,
           const Rcpp::NumericVector &curh0,      Rcpp::NumericMatrix &h0,
           const Rcpp::NumericMatrix &curpara,    Rcpp::NumericVector &para,
           const Rcpp::NumericVector &curlambda2, Rcpp::NumericMatrix &lambda2,
           const Rcpp::NumericMatrix &curtau2,    Rcpp::NumericVector &tau2,
           const Rcpp::NumericVector &curbeta, Rcpp::NumericVector &beta,
           const arma::umat &curmixind,  Rcpp::IntegerVector &mixind,
           const bool auxstore, const int thintime, const int where);

#endif
