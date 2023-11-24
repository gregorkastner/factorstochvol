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

#include <RcppArmadillo.h>
#include <stochvol.h>
#include "sampler.h"

using namespace Rcpp;
using namespace arma;

/*
 * As few as possible function arguments (better readable), or
 * many function arguments (than more has to be defined/initialized within)? 
 *  
 */

void update_fsv(arma::mat& facload,
                arma::mat& fac,
                arma::mat& logvar, //!!! if a factor is assumed to be homoskedastic (i.e. the corresponding element in heteroskedastic is 'false'), the corresponding column in logvar has to be initialized with 0s, otherwise the sampler will crash!
                arma::vec& logvar0,
                arma::mat& svpara,
                arma::mat& tau2,
                arma::vec& lambda2, // !!! if(ngprior==true) {if (columnwise==true) of length 'r', else of length 'm'} else of length 0 !!! not necessarily needed (as initial value) for sampling, however we want to store posterior draws
                arma::umat& curmixind, // not necessarily needed (as initial value) for sampling, however we want to store posterior draws
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
                const int& interweaving, // !!! cf documentation, which interweaving strategy can be employed in your case (if all(heteroskedastic==true) best option is 4 (deep interweaving), if some factors are homoscedastic 4 not possible, hence use 3 (shallow interweaving))
                const stochvol::ExpertSpec_FastSV& expert_idi,
                const stochvol::ExpertSpec_FastSV& expert_fac,
                const std::vector<stochvol::PriorSpec>& prior_specs,
                const double& B011inv,
                const bool& samplefac,
                const bool& signswitch,
                const int& i// rep
){
  
  const int m = y.n_rows;
  const int r = facload.n_cols;
  const int mpr = m+r;
  const int T = y.n_cols;
  
  // "linearized residuals"
  // NOTE: "log", "square" are component-wise functions, '*' denotes matrix multiplication
  arma::mat armaynorm(size(y));
  if (r > 0) {
    armaynorm = log(square(y - facload * fac));
  } else {
    armaynorm = log(square(y) + offset);
  }
  arma::mat armafnorm = log(square(fac));
  
  arma::mat armahtilde = exp(-logvar(arma::span::all, arma::span(0,m-1))/2.);
  
  
  // STEP 1:
  // update indicators, latent volatilities, and SV-parameters
  
  // STEP 1 for "linearized residuals"
  for (int j = 0; j < m; j++) {
    if (heteroskedastic(j) == true) {
      double curh0j = logvar0(j);
      arma::vec curh_j = logvar.unsafe_col(j);
      arma::uvec curmixind_j = curmixind.unsafe_col(j);
      double mu = svpara.at(0, j),
        phi = svpara.at(1, j),
        sigma = svpara.at(2, j);
      stochvol::update_fast_sv(armaynorm.row(j).t(), mu, phi, sigma, curh0j, curh_j, curmixind_j, prior_specs[j], expert_idi);
      svpara.at(0, j) = mu;
      svpara.at(1, j) = phi;
      svpara.at(2, j) = sigma;
      logvar0(j) = curh0j;
    } else {
      double rss;
      if (r > 0) {
        rss = sum(square(y.row(j) - facload.row(j)*fac));
      } else {
        rss = sum(square(y.row(j)));
      }
      const double sigma = 1/R::rgamma(priorhomoskedastic(j, 0) + .5*T, 1/(priorhomoskedastic(j, 1) + .5*rss));
      logvar.col(j).fill(log(sigma));
    }
  }
  
  // STEP 1 for factors
  for (int j = m; j < mpr; j++) {
    if (heteroskedastic(j) == true) {
      double curh0j = logvar0(j);
      arma::vec curh_j = logvar.unsafe_col(j);
      arma::uvec curmixind_j = curmixind.unsafe_col(j);
      double mu = 0,  //svpara.at(0, j),
        phi = svpara.at(1, j),
        sigma = svpara.at(2, j);
      stochvol::update_fast_sv(armafnorm.row(j-m).t(), mu, phi, sigma, curh0j, curh_j, curmixind_j, prior_specs[j], expert_fac);
      svpara.at(0, j) = 0;
      svpara.at(1, j) = phi;
      svpara.at(2, j) = sigma;
      logvar0(j) = curh0j;
    }
  }
  
  // intermediate step: calculate transformation of curh
  armahtilde = exp(-logvar(arma::span::all, arma::span(0,m-1))/2.);
  
  // STEP 2:
  // update factor loadings: m independent r-variate regressions
  // with T observations (for unrestricted case)
  if (r > 0) {
    // NEW: shrinkage part:
    // should we employ the NG-prior?
    if (ngprior) {
      if (!columnwise) {
        for (int j = 0; j < m; j++) {
          
          // draw lambda^2
          lambda2(j) = as<double>(rgamma(1, cShrink[j] + aShrink[j] * nonzerosperrow[j],
                                      1./(dShrink[j] + (aShrink[j]/2.) * sum(tau2.row(j)))));
          
          // draw tau^2
          for (int k = 0; k < r; k++) {
            if (restriction(j,k) != 0) {
              tau2(j,k) = do_rgig1(aShrink(j) - .5, facload(j,k)*facload(j,k),
                       aShrink(j)*lambda2(j));
            }
          }
        }
      } else {  // columnwise shrinkage
        for (int j = 0; j < r; j++) {
          
          // draw lambda^2
          lambda2(j) = as<double>(rgamma(1, cShrink[j] + aShrink[j] * nonzerospercol[j],
                                      1./(dShrink[j] + (aShrink[j]/2.) * sum(tau2.col(j)))));
          
          // draw tau^2
          for (int k = 0; k < m; k++) {
            if (restriction(k,j) != 0) {
              tau2(k,j) = do_rgig1(aShrink(j) - .5, facload(k,j)*facload(k,j),
                       aShrink(j)*lambda2(j));
            }
          }
        }
      }
    }
    
    /*
    // should we employ the DL-prior?
     
    //incorrect algorithm (wrong order of updating the hyperparameters)
    int unrestrictedelementcount = arma::accu(restriction); // would have to be initialized outside loop
    double tmpcounter4samplingtauDL; // would have to be initialized outside loop
    arma::mat armapsiDL(m,r); 
    arma::mat armaTDL(m,r); 
    arma::mat armaphiDL(m,r); // would be needed as function argument
    
    double tauDL; // DL prior, 1. incorrect original algorithm, 2. tauDL must be function argument...
    if (dlprior) {
      tmpcounter4samplingtauDL = 0;
      for (int j = 0; j < r; j++) {
        for (int k = 0; k < m; k++) {
          if (restriction(k,j) != 0) {
            armapsiDL(k,j) = 1. / do_rgig1(-.5, 1, (facload(k,j)*facload(k,j)) / (tauDL * tauDL * armaphiDL(k,j) * armaphiDL(k,j)));
            tmpcounter4samplingtauDL += fabs(facload(k,j))/armaphiDL(k,j);
            armaTDL(k,j) = do_rgig1(aShrink[0] - 1., 2*fabs(facload(k,j)), 1);
          }
        }
      }
      //Rprintf("%f\n", tmpcounter4samplingtauDL);
      //if (tmpcounter4samplingtauDL < 1.) Rprintf("THIS: %f\n", armaphiDL[0,0]);
      tauDL = do_rgig1((aShrink[0] - 1.) * unrestrictedelementcount, 2. * tmpcounter4samplingtauDL, 1);
      armaphiDL = armaTDL / accu(armaTDL);
      tau2 = armapsiDL % armaphiDL % armaphiDL * tauDL * tauDL;
    }
    */
    
    arma::mat armaXt(r,T);
    arma::vec armafacloadtmp = arma::zeros<arma::vec>(facloadtunrestrictedelements.size());
    arma::mat armaSigma(r,r);
    arma::mat armaR(r,r);
    arma::mat armaRinv(size(armaR));
    arma::mat armafacloadt = arma::trans(facload);
    int oldpos = 0;
    
    for (int j = 0; j < m; j++) {
      
      // TODO: some things outside
      
      // transposed design matrix Xt is filled "manually"
      int activecols = 0;
      for (int l = 0; l < r; l++) {
        if (restriction(j, l) != 0) {
          for (int k = 0; k < T; k++) {
            armaXt(activecols, k) = fac(l, k) * armahtilde(k, j);
          }
          activecols++;
        }
      }
      
      arma::colvec armaytilde = y.row(j).t() % armahtilde.col(j);
      
      
      // Now draw from the multivariate normal distribution
      // armaSigma is first used as temporary variable:
      armaSigma.submat(0,0,activecols-1,activecols-1) = armaXt.rows(0,activecols-1) * armaXt.rows(0,activecols-1).t();
      
      // add precisions to diagonal:
      armaSigma.submat(0,0,activecols-1,activecols-1).diag() += 1/arma::nonzeros(tau2.row(j));
      
      // Find Cholesky factor of posterior precision matrix
      try {
        armaR.submat(0, 0, activecols-1, activecols-1) = arma::chol(armaSigma.submat(0,0,activecols-1,activecols-1));
      } catch (...) {
        ::Rf_error("Error in run %i: Couldn't Cholesky-decompose posterior loadings precision in row %i", i+1, j+1);
      }
      
      
      // TODO: Check whether Armadillo automatically exploits the fact that R2 is upper triangular for inversion
      // (Partial) Answer: Seems to be OK for native R but solve(trimatu(R), I) is faster with OpenBLAS
      try {
        // armaRinv.submat(0,0,activecols-1,activecols-1) = arma::inv(arma::trimatu(armaR.submat(0,0,activecols-1,activecols-1)));
        armaRinv.submat(0,0,activecols-1,activecols-1) =
          arma::solve(arma::trimatu(armaR.submat(0,0,activecols-1,activecols-1)),
                      arma::eye<arma::mat>(activecols, activecols));
      } catch (...) {
        ::Rf_error("Error in run %i: Couldn't invert Cholesky factor of posterior loadings precision in row %i", i+1, j+1);
      }
      
      // calculate posterior covariance armaSigma:
      armaSigma.submat(0, 0, activecols-1, activecols-1) =
        armaRinv.submat(0, 0, activecols-1, activecols-1) *
        armaRinv.submat(0, 0, activecols-1, activecols-1).t();
      
      // calculate posterior mean:
      arma::colvec armamean(r);
      armamean.head(activecols) = armaSigma.submat(0, 0, activecols-1, activecols-1) *
        armaXt.submat(0, 0, activecols-1, T-1) *
        armaytilde;
      
      // draw from the r-variate normal distribution
      
      arma::colvec armadraw = rnorm(r);
      
      try {
        armafacloadtmp(arma::span(oldpos, oldpos + activecols - 1)) = armamean.head(activecols) + armaRinv.submat(0,0,activecols-1,activecols-1) * armadraw.head(activecols);
      } catch(...) {
        ::Rf_error("Error in run %i: Couldn't sample row %i of factor loadings", i+1, j+1);
      }
      
      //  Rprintf("\n%i to %i: ", oldpos, oldpos+activecols-1);
      //for (int is = oldpos; is < oldpos+activecols; is++) Rprintf("%f ", armafacloadtmp(is));
      //Rprintf("\n\n");
      oldpos = oldpos + activecols;
      
    }
    armafacloadt(facloadtunrestrictedelements) = armafacloadtmp;
    facload = arma::trans(armafacloadt);
    
    if (facloadtol > 0) {
      
      for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < r; jj++) {
          if (restriction(ii, jj) != 0) {
            if (facload(ii, jj) < facloadtol && facload(ii, jj) > 0.) {
              facload(ii, jj) = facloadtol;
            } else if (facload(ii, jj) > -facloadtol && facload(ii, jj) < 0.) {
              facload(ii, jj) = -facloadtol;
            } else if (facload(ii, jj) == 0.) {
              if (R::rbinom(1, 0.5) == 0) {
                facload(ii, jj) = facloadtol;
              } else {
                facload(ii, jj) = -facloadtol;
              }
            }
          }
        }
      }
    }
    
    //Rprintf("\n\n");
    //for (int is = 0; is < m; is++) Rprintf("%f %f\n", curfacload(is, 0), curfacload(is, 1));
    
    if (interweaving == 1 || interweaving == 3 || interweaving == 5 || interweaving == 7) {
      // STEP 2*: "Shallow" Interweaving
      
      //   // intermediate step: calculate transformation of curh
      //   armahtilde2 = exp(-logvar(arma::span::all, arma::span(m, m+r-1)));
      
      
      for (int j = 0; j < r; j++) {
        
        int userow = j;
        if (interweaving == 3 || interweaving == 7) { // find largest absolute element in column to interweave
          userow = 0;
          for (int k = 1; k < m; k++) if (fabs(facload(k, j)) > fabs(facload(userow, j))) userow = k;
        }
        if (interweaving == 5) { // find random nonzero element in column to interweave
          for (int k = 1; k < m; k++) {
            userow = floor(R::runif(0, m));
            if (fabs(facload(userow, j)) > 0.01) break;
          }
        }
        
        
        double newdiag2 = do_rgig1((nonzerospercol(j)- T) / 2.,
                                   sum(square(fac.row(j) * facload(userow,j))/exp(logvar.col(m+j)).t()),
                                   1/tau2(userow,j) + sum(square(nonzeros(facload.col(j))) / nonzeros(tau2.col(j))) / (facload(userow,j) * facload(userow,j)));
        
        double tmp = sqrt(newdiag2)/facload(userow,j);
        
        facload.col(j) *= tmp;
        fac.row(j) *= 1/tmp;
      }
    }
    
    if (interweaving == 2 || interweaving == 4 || interweaving == 6 || interweaving == 7) { // STEP 2+: "Deep" Interweaving
      for (int j = 0; j < r; j++) {
        
        int userow = j;
        if (interweaving == 4 || interweaving == 7) { // find largest absolute element in column to interweave
          userow = 0;
          for (int k = 1; k < m; k++) if (fabs(facload(k, j)) > fabs(facload(userow, j))) userow = k;
        }
        if (interweaving == 6) { // find random nonzero element in column to interweave
          for (int k = 1; k < m; k++) {
            userow = floor(R::runif(0, m));
            if (fabs(facload(userow, j)) > 0.01) break;
          }
          //Rprintf("use: %i ", userow);
        }
        
        
        //Rprintf("%i and %i\n", j, userow);
        
        double phi = svpara(1,m+j);//curpara(1,m+j);
        double sigma = svpara(2,m+j);
        double mu_old = log(facload(userow,j) * facload(userow,j));
        arma::vec hopen = logvar.col(m+j) + mu_old;
        double h0open = logvar0(m+j) + mu_old;
        double logacceptrate;
        double mu_prop;
        
        if (priorh0(m+j) < 0.) {  // old prior for h0 (stationary distribution, depends on phi), as in JCGS submission Feb 2016
          double tmph = hopen(0) - phi*h0open;
          for (int k = 1; k < T; k++) tmph += hopen(k) - phi*hopen(k-1);
          
          double gamma_old = (1 - phi) * mu_old;
          double gamma_prop = as<double>(rnorm(1, tmph/(T+B011inv), sigma/sqrt(T+B011inv)));
          mu_prop = gamma_prop/(1-phi);
          
          logacceptrate = logdnormquot(mu_prop, mu_old, h0open, sigma/sqrt(1-phi*phi));
          logacceptrate += logspecialquot(gamma_prop, gamma_old, .5, 1/(2.*tau2(userow,j)), 1-phi);
          logacceptrate += logdnormquot(gamma_old, gamma_prop, 0., sigma*sqrt(1/B011inv));
          
        } else {  // new prior does not depend on phi
          double tmph = hopen(0);
          for (int k = 1; k < (T-1); k++) tmph += hopen(k);
          
          double tmp4prop = T*priorh0(m+j)*(1-phi)*(1-phi) + 1;
          double prop_mean = (priorh0(m+j) * (1-phi) * (hopen(T-1) + (1-phi)*tmph - phi*h0open) + h0open) / tmp4prop;
          double prop_sd = (sqrt(priorh0(m+j)) * sigma) / sqrt(tmp4prop);
          
          mu_prop = as<double>(rnorm(1, prop_mean, prop_sd));
          logacceptrate = .5 * ((mu_prop - mu_old) - (exp(mu_prop) - exp(mu_old)) / tau2(userow,j));
        }
        
        // NEW, same for both priors:
        arma::vec relevantload = facload.col(j);
        arma::vec relevanttau2 = tau2.col(j);
        
        // use all except interwoven element (restricted loadings are assumed to be zero!)
        double mysum = accu(square(nonzeros(relevantload))/nonzeros(relevanttau2)) -
          (relevantload(userow)*relevantload(userow))/relevanttau2(userow);
        
        logacceptrate += .5 * ((nonzerospercol(j)-1)*(mu_prop - mu_old) -
          mysum / (facload(userow,j)*facload(userow,j)) * (exp(mu_prop) - exp(mu_old)));
        
        
        // Rprintf("ACCEPT? ");
        
        //ACCEPT/REJECT
        if (log(::unif_rand()) < logacceptrate) {
          //    Rprintf("ACC col %i el %02i - ", j+1, userow+1);
          //curh(_, m+j) = hopen - mu_prop;
          logvar.col(m+j) = hopen - mu_prop;
          logvar0(m+j) = h0open - mu_prop;
          
          double tmp = exp(mu_prop/2)/facload(userow,j);
          facload.col(j) *= tmp;
          fac.row(j) *= 1/tmp;
          //    } else {
          //     Rprintf("REJ col %i el %02i - ", j+1, userow+1);
        }
      }
      //   Rprintf("\n");
    }
    
    
    // STEP 3:
    // update the factors (T independent r-variate regressions with m observations)
    
    if (samplefac) {
      arma::colvec armadraw2(r*T);
      armadraw2.imbue(R::norm_rand); // .imbue: fill with values provided by a functor or lambda function; R::norm_rand is underlying function where R::rnorm() and Rcpp::rnorm() are wrappers around;
      arma::mat armaXt2(r,m);
      arma::mat armaSigma2(r,r);
      arma::mat armaR2(r,r);
      arma::mat armaR2inv(r,r);
      
      for (int j = 0; j < T; j++) {
        
        // transposed design matrix Xt2 (r x m) is filled "manually"
        for (int k = 0; k < m; k++) {
          for (int l = 0; l < r; l++) {
            armaXt2(l, k) = facload(k, l) * armahtilde(j,k);
          }
        }
        
        arma::colvec armaytilde2 = y.col(j) % armahtilde.row(j).t();
        
        // Now draw form the multivariate normal distribution
        
        // armaSigma2 is first used as temporary variable (to hold the precision):
        armaSigma2 = armaXt2 * armaXt2.t();
        
        // add precisions to diagonal:
        armaSigma2.diag() += exp(-logvar(j, arma::span(m, mpr-1)));
        
        // find Cholesky factor of posterior precision
        try {
          armaR2 = arma::chol(armaSigma2);
        } catch (...) {
          ::Rf_error("Error in run %i: Couldn't Cholesky-decompose posterior factor precision at time %i of %i", i+1, j+1, T);
        }
        
        try {
          //   armaR2inv = arma::inv(R2); # This is a little bit faster for very small matrices but a lot slower for large ones...
          //   armaR2inv = arma::inv(arma::trimatu(armaR2)); # This is OK on Native R but not so nice in OpenBLAS
          armaR2inv = arma::solve(arma::trimatu(armaR2), arma::eye<arma::mat>(r, r));
        } catch (...) {
          ::Rf_error("Error in run %i: Couldn't invert Cholesky factor of posterior factor precision at time %i of %i", i+1, j+1, T);
        }
        
        // calculate posterior covariance matrix armaSigma2:
        armaSigma2 = armaR2inv * armaR2inv.t();
        
        // calculate posterior mean armamean2:
        arma::colvec armamean2 = armaSigma2 * armaXt2 * armaytilde2;
        
        // draw from the r-variate normal distribution
        try {
          fac.col(j) = armamean2 + (armaR2inv * armadraw2.subvec(j*r, (j+1)*r - 1));
        } catch(...) {
          ::Rf_error("Error in run %i: Couldn't sample factors at time %i of %i", i+1, j+1, T);
        }
      }
    }
     
  }
  
  
  // SIGN SWITCH:
  if (signswitch) {
    for (int j = 0; j < r; j++) {
      if (as<double>(runif(1)) > .5) {
        facload.col(j) *= -1;
        fac.row(j) *= -1;
      }
    }
  }
  
}
