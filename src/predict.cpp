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

#include "predict.h"

using namespace Rcpp;

RcppExport SEXP predict(const SEXP obj_in, const SEXP store_in, const SEXP each_in) {

 // note: SEXP to Rcpp conversion REUSES memory unless "clone"d
 // Rcpp to Armadillo conversion allocates NEW memory unless deact'd

 const List obj(obj_in);

 const IntegerVector store(store_in);
 const int storelen = store.size();
 const int horizon = store(storelen - 1);
 const int each = as<int>(each_in);

 NumericVector h = obj["logvar"];
 const IntegerVector hDims = h.attr("dim");

 NumericVector facload = obj["facload"];
 const IntegerVector facloadDims = facload.attr("dim");

 const int m = facloadDims(0);
 const int r = facloadDims(1);
 const int mpr = m + r;
 const int len = facloadDims(2);

 NumericVector para = obj["para"];
 const IntegerVector paraDims = para.attr("dim");
 const arma::cube paras(para.begin(), paraDims(0), paraDims(1), paraDims(2), false);

 const arma::mat mus = paras.subcube(arma::span(0), arma::span(), arma::span());
 const arma::mat phis = paras.subcube(arma::span(1), arma::span(), arma::span());
 const arma::mat sigmas = paras.subcube(arma::span(2), arma::span(), arma::span());

 const arma::cube facloads(facload.begin(), facloadDims(0), facloadDims(1), facloadDims(2), false);

 const arma::cube hs(h.begin(), hDims(0), hDims(1), hDims(2), false);

 NumericVector hpreds_(mpr * len * each, 0.0);
 arma::cube hpredscub(hpreds_.begin(), mpr, len, each, false);
 arma::mat hpredsmat(hpreds_.begin(), mpr, len * each, false);

 for (int i = 0; i < each; i++) {
  hpredscub.slice(i) = hs.subcube(arma::span(hDims(0) - 1), arma::span(), arma::span());
 }

 NumericVector facpreds_(r * len * each, 0.0);
 arma::mat facpreds(facpreds_.begin(), r, len * each, false);

 NumericVector meanpreds_(m * len * each, 0.0);
 arma::mat meanpreds(meanpreds_.begin(), m, len * each, false);

 NumericVector volpredsstore_(m * len * each * storelen, 0.0);
 arma::cube volpredsstore(volpredsstore_.begin(), m, len * each, storelen, false);

 NumericVector meanpredsstore_(m * len * each * storelen, 0.0);
 arma::cube meanpredsstore(meanpredsstore_.begin(), m, len * each, storelen, false);

 arma::mat X(mpr, len);
 arma::mat Y(r, len * each);

 int storeindex = 0;

 //RNGScope scope;
 GetRNGstate(); // "by hand" because RNGScope isn't safe if return
                // variables are declared afterwards


 for (int i = 0; i < horizon; i++) {
  for (int e = 0; e < each; e++) {
   std::generate(X.begin(), X.end(), ::norm_rand); //R::norm_rand is underlying function where R::rnorm() and Rcpp::rnorm() are wrappers around;
   hpredscub.slice(e) = mus + phis % (hpredscub.slice(e) - mus) + sigmas % X;
  }
  if (r > 0) {
   std::generate(Y.begin(), Y.end(), ::norm_rand);
   facpreds = exp(hpredsmat.rows(m, m+r-1)/2) % Y;
   for (int j = 0; j < len * each; j++) {
    meanpreds.col(j) = facloads.slice(j/each) * facpreds.col(j);
   }
  }

  if (i == store(storeindex)-1) {
   volpredsstore.slice(storeindex) = exp(hpredsmat.rows(0, m-1)/2);
   if (r > 0) meanpredsstore.slice(storeindex) = meanpreds;
   storeindex++;
  }
 }

 List retval = List::create(
  Named("means") = wrap(meanpredsstore),
  Named("vols")  = wrap(volpredsstore)
 );

 PutRNGstate();
 return retval;
}
