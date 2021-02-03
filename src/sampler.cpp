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

#include <RcppArmadillo.h>
#include "sampler.h"

using namespace Rcpp;

// rgig is imported from GIGrvg
double rgig1(double lambda, double chi, double psi) {
 PutRNGstate();
 SEXP (*fun)(SEXP, SEXP, SEXP, SEXP) = NULL;
 if (!fun) {
  fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("GIGrvg", "rgig");
 }
 GetRNGstate();
 return as<double>(fun(wrap(1), wrap(lambda), wrap(chi), wrap(psi)));
}

double do_rgig1(double lambda, double chi, double psi) { 
 SEXP (*fun)(int, double, double, double) = NULL;
 if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
 return as<double>(fun(1, lambda, chi, psi));
}

void test(double * data, int size) {
 for (int i = 0; i < size; i++) {
  data[i] = i;
 }
}

RcppExport SEXP sampler(const SEXP y_in, const SEXP draws_in,
  const SEXP burnin_in, const SEXP startval_in,
  const SEXP bmu_in, const SEXP Bmu_in,
  const SEXP priorphi_in, const SEXP Bsigma_in,
  const SEXP priorbeta_in,
  const SEXP model_mean_in, const SEXP shrinkagepriors_in,
  const SEXP thin_in, const SEXP auxstore_in,
  const SEXP thintime_in,
  const SEXP quiet_in, const SEXP para_in,
  const SEXP MHsteps_in, const SEXP B011_in, const SEXP B022_in,
  const SEXP mhcontrol_in, const SEXP gammaprior_in,
  const SEXP offset_in, const SEXP truncnormal_in,
  const SEXP restriction_in, const SEXP interweaving_in, 
  const SEXP signswitch_in, const SEXP runningstore_in,
  const SEXP runningstoreevery_in, const SEXP runningstoremoments_in,
  const SEXP pfl_in, const SEXP heteroskedastic_in,
  const SEXP priorhomoskedastic_in, const SEXP priorh0_in,
  const SEXP samplefac_in) {

 // note: SEXP to Rcpp conversion REUSES memory unless "clone"d
 // Rcpp to Armadillo conversion allocates NEW memory unless deact'd
 NumericMatrix y(y_in);
 List startval(clone(startval_in));
 
 // interweaving strategy:
 // 0 = none
 // 1 = "shallow interweaving" (diagonal element)
 // 2 = "deep interweaving" (diagonal element)
 // 3 = "shallow interweaving" (largest |element|)
 // 4 = "deep interweaving" (largest |element|)
 // 5 = "shallow interweaving" (random element)
 // 6 = "deep interweaving" (random element)
 
 const int interweaving = as<int>(interweaving_in);

 const bool signswitch  = as<bool>(signswitch_in);
 const int runningstore = as<int>(runningstore_in);
 
 const bool samplefac   = as<bool>(samplefac_in);

 const bool model_mean  = as<bool>(model_mean_in);

 /*
  * LOOP STORAGE (current_* variables)
  */

 //current factor loadings matrix draws
 NumericMatrix curfacload = startval["facload"];
 
 const int m = y.nrow(); // number of time series
 const int T = y.ncol(); // length of time series
 const int r = curfacload.ncol(); // number of latent factors
 const int mpr = m + r;

 // 1 = Normal, 2 = NG (rowwise), 3 = NG (colwise), 4 = DL
 const int pfl = as<int>(pfl_in);
 bool ngprior = false;
 if (r > 0 && (pfl == 2 || pfl == 3)) ngprior = true;
 bool columnwise = false;
 if (r > 0 && pfl == 3) columnwise = true;
 bool dlprior = false;
 if (r > 0 && pfl == 4) dlprior = true;
 NumericVector sv(heteroskedastic_in);
 NumericMatrix priorhomoskedastic(priorhomoskedastic_in);
 NumericVector priorh0(priorh0_in);
 IntegerMatrix restr(restriction_in);
 arma::imat armarestr(restr.begin(), restr.nrow(), restr.ncol(), false);
 
 arma::irowvec nonzerospercol = arma::sum(armarestr, 0);
 arma::icolvec nonzerosperrow = arma::sum(armarestr, 1);

 // restriction on factor loadings matrix:
 for (int i = 0; i < curfacload.nrow(); i++) {
  for (int j = 0; j < curfacload.ncol(); j++) {
   if (armarestr(i, j) == 0) curfacload(i,j) = 0.;
  }
 }

 //convention: "arma"-prefixed variables denote Armadillo proxy objects
 const arma::mat armay_original(y.begin(), m, T, false);
 arma::mat armay = armay_original;  //(armay_original.begin(), m, T, model_mean);  // demeaned
 arma::mat armay_regression = armay;
 arma::mat armafacload(curfacload.begin(), m, r, false);
 arma::mat armafacloadt = arma::trans(armafacload);

 arma::uvec armafacloadtunrestrictedelements = arma::find(armarestr.t() != 0);
 //for (int i = 0; i < 10; i++) Rprintf("%i ", armafacloadtunrestrictedelements(i));
 //Rprintf("\n\n");
 arma::vec armafacloadtmp = arma::zeros<arma::vec>(armafacloadtunrestrictedelements.size());

 // curfacloadinter will hold the diagonal factor loadings stemming from the interwoven sampler
 NumericVector curfacload2inter(r);
 arma::vec armafacload2inter(curfacload2inter.begin(), curfacload2inter.length(), false);

 //current factor draws
 NumericMatrix curf = startval["fac"];
 arma::mat armaf(curf.begin(), curf.nrow(), curf.ncol(), false);
 
// //current "centered" factor loadings
// NumericMatrix curfstar(curf.nrow(), curf.ncol());
// arma::mat armafstar(curfstar.begin(), curfstar.nrow(), curfstar.ncol(), false);

 //current log-volatility draws
 NumericMatrix curh = startval["latent"]; // does not contain h0!
 arma::mat armah(curh.begin(), curh.nrow(), curh.ncol(), false);

 //transformation thereof (cols 1 to m)
 NumericMatrix curhtilde(curh.nrow(), m);

 arma::mat armahtilde(curhtilde.begin(), curhtilde.nrow(), curhtilde.ncol(), false);

// //transformation thereof (cols m+1 to m+r)
// NumericMatrix curhtilde2(curh.nrow(), r);
// arma::mat armahtilde2(curhtilde2.begin(), curhtilde2.nrow(), curhtilde2.ncol(), false);
 
 NumericVector curh0 = startval["latent0"];
 arma::vec armah0(curh0.begin(), curh0.length(), false);
 
 //current parameter draws
 const List startpara = startval["para"];
 const NumericVector startmu = startpara["mu"];
 const NumericVector startphi = startpara["phi"];
 const NumericVector startsigma = startpara["sigma"];
 
 NumericMatrix curpara(3, mpr);
 for (int i = 0; i < m; i++) {
  curpara(0,i) = startmu(i);
  curpara(1,i) = startphi(i);
  curpara(2,i) = startsigma(i);
 }

 for (int i = m; i < mpr; i++) {
  curpara(0,i) = 0.;
  curpara(1,i) = startphi(i);
  curpara(2,i) = startsigma(i);
 }

 //current mixture indicator draws
 arma::umat curmixind(T, mpr);
 
 //current mixture probability draws
 NumericVector curmixprob(10 * T * mpr);
 
 // shrinkage prior:
 const List shrinkagepriors(shrinkagepriors_in);

 // NA means: use N(0,tau2)-prior with tau2 fixed
 const NumericVector aShrink = shrinkagepriors["a"];
 const NumericVector cShrink = shrinkagepriors["c"];
 const NumericVector dShrink = shrinkagepriors["d"];


 int nlambda;
 if (ngprior) {
  if (columnwise) {
   nlambda = r;
  } else {
   nlambda = m;
  }
 } else nlambda = 0;

 //current shrinkage latents lambda^2
 NumericVector curlambda2(nlambda);
 arma::vec armalambda2(curlambda2.begin(), curlambda2.size(), false);
 
 //current shrinkage variances tau^2
 NumericMatrix curtau2 = startval["tau2"];

 // current regression betas
 // note that only single regression is implemented
 NumericVector curbeta(m); curbeta.fill(0);
 arma::vec armabeta(curbeta.begin(), curbeta.size(), false);

 /*
  * MARKOV CHAIN AND MODEL SETUP
  */
 
 // restriction on factor loadings matrix:
 for (int i = 0; i < curtau2.nrow(); i++) {
  for (int j = 0; j < curtau2.ncol(); j++) {
   if (armarestr(i,j) == 0) curtau2(i,j) = 0.;
  }
 }
 
 int unrestrictedelementcount = arma::accu(armarestr);
 arma::mat armatau2(curtau2.begin(), curtau2.nrow(), curtau2.ncol(), false);
 double tauDL = 1.; // whatever
 double tmpcounter4samplingtauDL;
 arma::mat armapsiDL(m, r, arma::fill::zeros);  // used for DL prior
 arma::mat armaphiDL(m, r, arma::fill::zeros);  // used for DL prior
 arma::mat armaTDL(m, r, arma::fill::zeros);  // used for DL prior
 for (int i = 0; i < m; i++) {
  for (int j = 0; j < r; j++) {
   if (armarestr(i,j) != 0) {
    armapsiDL(i,j) = 1./unrestrictedelementcount;
    armaphiDL(i,j) = 1./unrestrictedelementcount;
    armaTDL(i,j) = 1./unrestrictedelementcount;
   }
  }
 }

 // number of MCMC draws
 const int burnin = as<int>(burnin_in);
 const int draws  = as<int>(draws_in);
 const int N 	    = burnin + draws;
 
 // temporary stroage for hopen in interweaving
 NumericVector hopen(T);

 // prior parameters
 const double bmu    = as<double>(bmu_in);
 const double Bmu    = as<double>(Bmu_in);
 
 const NumericVector priorphi(priorphi_in);
 const double a0idi     = priorphi(0);
 const double b0idi     = priorphi(1);
 const double a0fac     = priorphi(2);
 const double b0fac     = priorphi(3);
 const NumericVector Bsigma(Bsigma_in);

 const NumericVector priorbeta(priorbeta_in);
 
 //std::fill(armatau2.begin(), armatau2.end(), 1.);
 
 // thinning parameters
 const int thintime  = as<int>(thintime_in);
 const bool auxstore = as<bool>(auxstore_in);
 const int thin      = as<int>(thin_in);
 const int runningstoreevery = as<int>(runningstoreevery_in);
 const int runningstoremoments = as<int>(runningstoremoments_in);

 // verbosity control
 const bool verbose = !as<bool>(quiet_in);

 // "expert" settings:
 const double B011inv         = 1./as<double>(B011_in);
 const double B022inv         = 1./as<double>(B022_in);
 const bool Gammaprior        = as<bool>(gammaprior_in);
 const bool truncnormal       = as<bool>(truncnormal_in);
 const double MHcontrol       = as<double>(mhcontrol_in);
 const int MHsteps            = as<int>(MHsteps_in);
 const int parameterization   = as<int>(para_in);
 const stochvol::ExpertSpec_FastSV expert_idi {
   parameterization > 2,  // interweave
   stochvol::Parameterization::CENTERED,  // centered_baseline always
   B011inv,
   B022inv,
   MHsteps,
   MHcontrol < 0 ? stochvol::ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE : stochvol::ExpertSpec_FastSV::ProposalSigma2::LOG_RANDOM_WALK,
   MHcontrol,
   truncnormal ? stochvol::ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL : stochvol::ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL
 };
 const stochvol::ExpertSpec_FastSV expert_fac {
   parameterization > 2,  // interweave
   stochvol::Parameterization::CENTERED,  // centered_baseline always
   B011inv,
   B022inv,
   3,
   MHcontrol < 0 ? stochvol::ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE : stochvol::ExpertSpec_FastSV::ProposalSigma2::LOG_RANDOM_WALK,
   MHcontrol,
   truncnormal ? stochvol::ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL : stochvol::ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL
 };

 // offset (only used in zero-factor model)
 const double offset          = as<double>(offset_in);

 // moment-matched IG-prior
 const NumericVector C0 = 1.5*Bsigma;

 // prior specification object for stochvol
 std::vector<stochvol::PriorSpec> prior_specs(mpr);
 {
   using stochvol::PriorSpec;
   for (int j = 0; j < m; j++) {
     prior_specs[j] = {
       (priorh0(j) <= 0) ? PriorSpec::Latent0() : PriorSpec::Latent0(PriorSpec::Constant(priorh0(j))),
         PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),
         PriorSpec::Phi(PriorSpec::Beta(a0idi, b0idi)),
         Gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma(j))) : PriorSpec::Sigma2(PriorSpec::InverseGamma(2.5, C0(j))),
         PriorSpec::Nu(PriorSpec::Infinity()),
         PriorSpec::Rho(PriorSpec::Constant(0)),
         PriorSpec::Covariates(PriorSpec::MultivariateNormal{{priorbeta[0]}, {std::pow(priorbeta[1], -2)}})
     };
   }
   for (int j = m; j < mpr; j++) {
     prior_specs[j] = {
       (priorh0(j) <= 0) ? PriorSpec::Latent0() : PriorSpec::Latent0(PriorSpec::Constant(priorh0(j))),
         PriorSpec::Mu(PriorSpec::Constant(0)),
         PriorSpec::Phi(PriorSpec::Beta(a0fac, b0fac)),
         Gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma(j))) : PriorSpec::Sigma2(PriorSpec::InverseGamma(2.5, C0(j)))
     };
   }
 }
 
 /* 
  * FINAL STORAGE (returned to R)
  */

 // NOTE: (Almost) all storage of MCMC draws is done in NumericVectors
 // because no 'array' structure is available at this point in time.
 
 // facload holds the factor loadings:
 NumericVector facload(curfacload.nrow() * curfacload.ncol() * (draws/thin));
 facload.attr("dim") = Dimension(curfacload.nrow(), curfacload.ncol(), draws/thin);

 
 int timestores = 0;
 if (thintime == -1) {
  timestores = 1;  // keep only last vols/facs
 } else if (thintime == 1) {  // keep all latent vols/facs
  timestores = T;
 } else if (thintime > 1) {  // keep some
  timestores = T/thintime;
 } 

 // h holds the latent log-volatilities, but not h0!
 NumericVector h(timestores * curh.ncol() * (draws/thin));
 h.attr("dim") = Dimension(timestores, curh.ncol(), draws/thin);

 // f holds the latent factor draws
 NumericVector f(curf.nrow() * timestores * (draws/thin));
 f.attr("dim") = Dimension(curf.nrow(), timestores, draws/thin);
 
 long tmplength;  // don't need to allocate much if not used
 if (runningstore >= 1) tmplength = T; else tmplength = 0;

 // hrunmean holds the running mean of the latent log-volatilities
 NumericMatrix hrunmean(tmplength, curh.ncol());
 arma::mat armahrunmean(hrunmean.begin(), hrunmean.nrow(), hrunmean.ncol(), false);
 
 // hrunm2 holds the running second moments of the latent log-volatilities
 NumericMatrix hrunm2(tmplength*(runningstoremoments >= 2), curh.ncol());
 arma::mat armahrunm2(hrunm2.begin(), hrunm2.nrow(), hrunm2.ncol(), false);
 
 // hrunm3 holds the running third moments of the latent log-volatilities
 NumericMatrix hrunm3(tmplength*(runningstoremoments >= 3), curh.ncol());
 arma::mat armahrunm3(hrunm3.begin(), hrunm3.nrow(), hrunm3.ncol(), false);
 
 // hrunm4 holds the running fourth moments of the latent log-volatilities
 NumericMatrix hrunm4(tmplength*(runningstoremoments >= 4), curh.ncol());
 arma::mat armahrunm4(hrunm4.begin(), hrunm4.nrow(), hrunm4.ncol(), false);

// NumericMatrix hrunmin(tmplength, curh.ncol());
// hrunmin.fill(100000.);
// arma::mat armahrunmin(hrunmin.begin(), hrunmin.nrow(), hrunmin.ncol(), false);
 
// NumericMatrix hrunmax(tmplength, curh.ncol());
// hrunmax.fill(-100000.);
// arma::mat armahrunmax(hrunmax.begin(), hrunmax.nrow(), hrunmax.ncol(), false);

 if (runningstore >= 6) tmplength = T; else tmplength = 0;

 arma::mat curcom(tmplength, m+1);
 
 // comrunmean holds the running mean of the "communalities"
 NumericMatrix comrunmean(tmplength, m+1);
 arma::mat armacomrunmean(comrunmean.begin(), comrunmean.nrow(), comrunmean.ncol(), false);
 
 // comrunm2 holds the running second moments of the "communalities"
 NumericMatrix comrunm2(tmplength*(runningstoremoments >= 2), m+1);
 arma::mat armacomrunm2(comrunm2.begin(), comrunm2.nrow(), comrunm2.ncol(), false);
 
 // comrunm3 holds the running third moments of the "communalities"
 NumericMatrix comrunm3(tmplength*(runningstoremoments >= 3), m+1);
 arma::mat armacomrunm3(comrunm3.begin(), comrunm3.nrow(), comrunm3.ncol(), false);
 
 // comrunm4 holds the running fourth moments of the "communalities"
 NumericMatrix comrunm4(tmplength*(runningstoremoments >= 4), m+1);
 arma::mat armacomrunm4(comrunm4.begin(), comrunm4.nrow(), comrunm4.ncol(), false);

 if (runningstore >= 5) tmplength = T; else tmplength = 0;
 
 arma::mat tmpcor(T, m*(m-1)/2);
 arma::mat tmpsds(T, m);
 
 int tmpcounter = 0;
 arma::uvec diagindices(m);

 for (int k = 0; k < m; k++) {
  for (int l = k; l < m; l++) {
   if (k == l) diagindices(k) = tmpcounter;
   tmpcounter++;
  }
 }
 
 // holds the running mean of correlation matrix
 // NOTE: Manual storage of lower-triangular portion (excluding diagonal), column major!
 NumericMatrix corrunmean(tmplength, (m*(m-1))/2);
 arma::mat armacorrunmean(corrunmean.begin(), corrunmean.nrow(), corrunmean.ncol(), false);

 // holds the running second moments of correlation matrix
 NumericMatrix corrunm2(tmplength*(runningstoremoments >= 2), (m*(m-1))/2);
 arma::mat armacorrunm2(corrunm2.begin(), corrunm2.nrow(), corrunm2.ncol(), false);
 
 // holds the running third moments of correlation matrix
 NumericMatrix corrunm3(tmplength*(runningstoremoments >= 3), (m*(m-1))/2);
 arma::mat armacorrunm3(corrunm3.begin(), corrunm3.nrow(), corrunm3.ncol(), false);
 
 // holds the running fourth moments of correlation matrix
 NumericMatrix corrunm4(tmplength*(runningstoremoments >= 4), (m*(m-1))/2);
 arma::mat armacorrunm4(corrunm4.begin(), corrunm4.nrow(), corrunm4.ncol(), false);

 if (runningstore >= 4) tmplength = T; else tmplength = 0;
 
 arma::mat tmpcov(T, m*(m+1)/2);
 arma::mat tmpvol(T, m);

 // holds the running mean of covariance matrix
 // NOTE: Manual storage of lower-triangular portion (including diagonal), column major!
 NumericMatrix covrunmean(tmplength, (m*(m+1))/2);
 arma::mat armacovrunmean(covrunmean.begin(), covrunmean.nrow(), covrunmean.ncol(), false);

 // holds the running second moments of covariance matrix
 NumericMatrix covrunm2(tmplength*(runningstoremoments >= 2), (m*(m+1))/2);
 arma::mat armacovrunm2(covrunm2.begin(), covrunm2.nrow(), covrunm2.ncol(), false);
 
 // holds the running third moments of covariance matrix
 NumericMatrix covrunm3(tmplength*(runningstoremoments >= 3), (m*(m+1))/2);
 arma::mat armacovrunm3(covrunm3.begin(), covrunm3.nrow(), covrunm3.ncol(), false);
 
 // holds the running fourth moments of covariance matrix
 NumericMatrix covrunm4(tmplength*(runningstoremoments >= 4), (m*(m+1))/2);
 arma::mat armacovrunm4(covrunm4.begin(), covrunm4.nrow(), covrunm4.ncol(), false);

 // holds the running mean of sqrt(diag(covariance matrix)) ("volatilities")
 NumericMatrix volrunmean(tmplength, m);
 arma::mat armavolrunmean(volrunmean.begin(), volrunmean.nrow(), volrunmean.ncol(), false);

 // note: second and fourth moments of sqrt(diag(covariance matrix)) are
 // already stored in covrunmean and covrunm2, but for convenience reasons
 // we just do it again (comparably cheap)
 
 // holds the running second moments of sqrt(diag(covariance matrix)) ("volatilities")
  NumericMatrix volrunm2(tmplength*(runningstoremoments >= 2), m);
 arma::mat armavolrunm2(volrunm2.begin(), volrunm2.nrow(), volrunm2.ncol(), false);
 
 // holds the running third moments of sqrt(diag(covariance matrix)) ("volatilities")
 NumericMatrix volrunm3(tmplength*(runningstoremoments >= 3), m);
 arma::mat armavolrunm3(volrunm3.begin(), volrunm3.nrow(), volrunm3.ncol(), false);
 
 // holds the running fourth moments of sqrt(diag(covariance matrix)) ("volatilities")
 NumericMatrix volrunm4(tmplength*(runningstoremoments >= 4), m);
 arma::mat armavolrunm4(volrunm4.begin(), volrunm4.nrow(), volrunm4.ncol(), false);
 
 if (runningstore >= 3) tmplength = T; else tmplength = 0;

 arma::mat htranstmp(tmplength, m+r);
 
 // hrunmeantrans holds the running mean of exp(latent log-volatilities/2)
 NumericMatrix hrunmeantrans(tmplength, curh.ncol());
 arma::mat armahrunmeantrans(hrunmeantrans.begin(), hrunmeantrans.nrow(), hrunmeantrans.ncol(), false);

 // hrunm2trans holds the running second moments of exp(latent log-volatilities/2)
 NumericMatrix hrunm2trans(tmplength*(runningstoremoments >= 2), curh.ncol());
 arma::mat armahrunm2trans(hrunm2trans.begin(), hrunm2trans.nrow(), hrunm2trans.ncol(), false);
 
 // hrunm3trans holds the running third moments of exp(latent log-volatilities/2)
 NumericMatrix hrunm3trans(tmplength*(runningstoremoments >= 3), curh.ncol());
 arma::mat armahrunm3trans(hrunm3trans.begin(), hrunm3trans.nrow(), hrunm3trans.ncol(), false);
 
 // hrunm4trans holds the running fourth moments of exp(latent log-volatilities/2)
 NumericMatrix hrunm4trans(tmplength*(runningstoremoments >= 4), curh.ncol());
 arma::mat armahrunm4trans(hrunm4trans.begin(), hrunm4trans.nrow(), hrunm4trans.ncol(), false);
   
 if (runningstore >= 2) tmplength = T; else tmplength = 0;

 // frunmean hold the running mean of the latent factors
 NumericMatrix frunmean(curf.nrow(), tmplength);
 arma::mat armafrunmean(frunmean.begin(), frunmean.nrow(), frunmean.ncol(), false);

 // frunm2 holds the running second moments of the latent log-volatilities
 NumericMatrix frunm2(curf.nrow(), tmplength*(runningstoremoments >= 2));
 arma::mat armafrunm2(frunm2.begin(), frunm2.nrow(), frunm2.ncol(), false);
 
 // frunm3 holds the running third moments of the latent log-volatilities
 NumericMatrix frunm3(curf.nrow(), tmplength*(runningstoremoments >= 3));
 arma::mat armafrunm3(frunm3.begin(), frunm3.nrow(), frunm3.ncol(), false);
 
 // frunm4 holds the running fourth moments of the latent log-volatilities
 NumericMatrix frunm4(curf.nrow(), tmplength*(runningstoremoments >= 4));
 arma::mat armafrunm4(frunm4.begin(), frunm4.nrow(), frunm4.ncol(), false);
 
// NumericMatrix frunmin(curf.nrow(), tmplength);
// frunmin.fill(100000.);
// arma::mat armafrunmin(frunmin.begin(), frunmin.nrow(), frunmin.ncol(), false);
 
// NumericMatrix frunmax(curf.nrow(), tmplength);
// frunmax.fill(-100000.);
// arma::mat armafrunmax(frunmax.begin(), frunmax.nrow(), frunmax.ncol(), false);

 int auxstoresize;
 if (auxstore) {
  auxstoresize = (draws/thin);
 } else {
  auxstoresize = 0;  // just a dummy object
 }
 
 // mixind holds the mixture indicators for the auxiliary mixture sampling
 IntegerVector mixind(T * mpr * auxstoresize);
 if (auxstore) mixind.attr("dim") = Dimension(T, mpr, draws/thin);
 
 // mixprob holds the mixture probabilities for the auxmix
 //NumericVector mixprob(10 * T * mpr * auxstoresize);
 //mixprob.attr("dim") = Dimension(10, T, mpr, draws/thin); no 4-dim possible?

 // lambda2 holds the latent lambda^2 draws
 NumericMatrix lambda2(curlambda2.length(), auxstoresize);
 
 // tau2 holds the latent variances
 NumericVector tau2(curtau2.nrow() * curtau2.ncol() * auxstoresize);
 tau2.attr("dim") = Dimension(curtau2.nrow(), curtau2.ncol(), auxstoresize);
 
 // h0 holds the initival latent log-volatilities
 NumericMatrix h0(curh0.length(), draws/thin);
 
 // para holds the parameter draws (mu, phi, sigma)
 NumericVector para(3 * mpr * (draws/thin));
 para.attr("dim") = Dimension(3, mpr, draws/thin) ;

 // beta holds the beta draws
 NumericVector beta(model_mean * m * (draws/thin));
 beta.attr("dim") = Dimension(model_mean * m, draws/thin);

 /*
  * TEMPORARY STORAGE
  */
  
 // curynorm will hold log((y - facload %*% f)^2) in STEP 1
 NumericMatrix curynorm(y.nrow(), y.ncol());
 arma::mat armaynorm(curynorm.begin(), curynorm.nrow(), curynorm.ncol(), false);
 
 // curynorm2 will hold log(f^2) in STEP 1
 NumericMatrix curfnorm(curf.nrow(), curf.ncol());
 arma::mat armafnorm(curfnorm.begin(), curfnorm.nrow(), curfnorm.ncol(), false);

 // Xt will hold the transposed design matrix in STEP 2
 NumericMatrix Xt(r, T);
 arma::mat armaXt(Xt.begin(), Xt.nrow(), Xt.ncol(), false);
 
 // Xt2 will hold the transposed design matrix in STEP 3
 NumericMatrix Xt2(r, m);
 arma::mat armaXt2(Xt2.begin(), Xt2.nrow(), Xt2.ncol(), false);

 // ytilde will hold the normalized observation vector in STEP 2
 NumericVector ytilde(T);
 arma::colvec armaytilde(ytilde.begin(), ytilde.length(), false);
 
 // ytilde2 will hold the normalized observation vector in STEP 3
 NumericVector ytilde2(m);
 arma::colvec armaytilde2(ytilde2.begin(), ytilde2.length(), false);

 // Sigma will hold the posterior variance-covariance matrix in STEP 2
 NumericMatrix Sigma(r, r);
 arma::mat armaSigma(Sigma.begin(), Sigma.nrow(), Sigma.ncol(), false);
 
 // R will hold the Cholesky factor of posterior precision matrix in STEP 2
 NumericMatrix R(r, r);
 arma::mat armaR(R.begin(), R.nrow(), R.ncol(), false);
 
 // Rinv will hold the inverse Cholesky factor of posterior precision matrix in STEP 3
 NumericMatrix Rinv(r, r);
 arma::mat armaRinv(Rinv.begin(), Rinv.nrow(), Rinv.ncol(), false);

 // Sigma2 will hold the posterior variance-covariance matrix in STEP 3
 NumericMatrix Sigma2(r, r);
 arma::mat armaSigma2(Sigma2.begin(), Sigma2.nrow(), Sigma2.ncol(), false);
 
 // R2 will hold the Cholesky factor of posterior precision matrix in STEP 3
 NumericMatrix R2(r, r);
 arma::mat armaR2(R2.begin(), R2.nrow(), R2.ncol(), false);
 
 // R2inv will hold the inverse Cholesky factor of posterior precision matrix in STEP 3
 NumericMatrix R2inv(r, r);
 arma::mat armaR2inv(R2inv.begin(), R2inv.nrow(), R2inv.ncol(), false);

 // mean will hold the posterior mean in STEP 2
 NumericVector mean(r);
 arma::colvec armamean(mean.begin(), mean.length(), false);
 
 // mean2 will hold the posterior mean in STEP 3
 NumericVector mean2(r);
 arma::colvec armamean2(mean2.begin(), mean2.length(), false);

 // draw will hold some random draws in STEP 2 
 
// NumericVector draw(step2draws);
 NumericVector draw(r);
 arma::colvec armadraw(draw.begin(), draw.length(), false);

 // draw2 will hold some random draws in STEP 3
 NumericVector draw2(r*T);
 arma::colvec armadraw2(draw2.begin(), draw2.length(), false);

 int effi = -burnin;
 int effirunningstore = 0;

 // temporary variables for the updated stochvol code
 arma::mat curpara_arma(curpara.begin(), curpara.nrow(), curpara.ncol(), false);
 arma::mat curh_arma(curh.begin(), curh.nrow(), curh.ncol(), false);
 arma::vec beta_j(1);

 RNGScope scope;
 //GetRNGstate(); // "by hand" because RNGScope isn't safe if return
                // variables are declared afterwards

 for (int j = m; j < mpr; j++) {
  if (sv(j) == false) {
   armah.col(j).fill(0.);
  }
 }

 int space4print = floor(log10(N + .1)) + 1;
 int doevery = ceil((2000.*N)/((r+1)*T*m));

 for (int i = 0; i < N; i++, effi++) {  // BEGIN main MCMC loop
  
  if (verbose && (i % doevery == 0)) {
   Rprintf("\r********* Iteration %*i of %*i (%3.0f%%) *********",
    space4print, i+1, space4print, N, 100.*(i+1)/N);
  }

  if (i % 20 == 0) {
    ::R_CheckUserInterrupt();
  }

  // "linearized residuals"
  // NOTE: "log", "square" are component-wise functions, '*' denotes matrix multiplication
  if (r > 0) {
   armaynorm = log(square(armay - armafacload * armaf));
  } else {
   armaynorm = log(square(armay) + offset);
  }
  armafnorm = log(square(armaf));
  
  armahtilde = exp(-armah(arma::span::all, arma::span(0,m-1))/2.);


  // STEP 1:
  // update indicators, latent volatilities, and SV-parameters
  


  // STEP 1 for "linearized residuals"
  for (int j = 0; j < m; j++) {
   if (sv(j) == true) {
     double curh0j = curh0(j);
     arma::vec curh_j = armah.unsafe_col(j);
     arma::uvec curmixind_j = curmixind.unsafe_col(j);
     double mu = curpara_arma.at(0, j),
            phi = curpara_arma.at(1, j),
            sigma = curpara_arma.at(2, j);
     stochvol::update_fast_sv(armaynorm.row(j).t(), mu, phi, sigma, curh0j, curh_j, curmixind_j, prior_specs[j], expert_idi);
     curpara_arma.at(0, j) = mu;
     curpara_arma.at(1, j) = phi;
     curpara_arma.at(2, j) = sigma;
     curh0(j) = curh0j;
   } else {
     double rss;
     if (r > 0) {
       rss = sum(square(armay.row(j) - armafacload.row(j)*armaf));
     } else {
       rss = sum(square(armay.row(j)));
     }
     const double sigma = 1/R::rgamma(priorhomoskedastic(j, 0) + .5*T, 1/(priorhomoskedastic(j, 1) + .5*rss));
     armah.col(j).fill(log(sigma));
   }
  }

  // STEP 1 for factors
  for (int j = m; j < mpr; j++) {
   if (sv(j) == true) {
     double curh0j = curh0(j);
     arma::vec curh_j = armah.unsafe_col(j);
     arma::uvec curmixind_j = curmixind.unsafe_col(j);
     double mu = 0,  //curpara_arma.at(0, j),
            phi = curpara_arma.at(1, j),
            sigma = curpara_arma.at(2, j);
     stochvol::update_fast_sv(armafnorm.row(j-m).t(), mu, phi, sigma, curh0j, curh_j, curmixind_j, prior_specs[j], expert_fac);
     curpara_arma.at(0, j) = 0;
     curpara_arma.at(1, j) = phi;
     curpara_arma.at(2, j) = sigma;
     curh0(j) = curh0j;
   }
  }
  
  // intermediate step: calculate transformation of curh
  armahtilde = exp(-armah(arma::span::all, arma::span(0,m-1))/2.);

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
      armalambda2(j) = as<double>(rgamma(1, cShrink[j] + aShrink[j] * nonzerosperrow[j],
         1./(dShrink[j] + (aShrink[j]/2.) * sum(armatau2.row(j)))));
   
      // draw tau^2
      for (int k = 0; k < r; k++) {
       if (armarestr(j,k) != 0) {
	armatau2(j,k) = do_rgig1(aShrink(j) - .5, armafacload(j,k)*armafacload(j,k),
	 aShrink(j)*armalambda2(j));
       }
      }
     }
    } else {  // columnwise shrinkage
     for (int j = 0; j < r; j++) {

      // draw lambda^2
      armalambda2(j) = as<double>(rgamma(1, cShrink[j] + aShrink[j] * nonzerospercol[j],
         1./(dShrink[j] + (aShrink[j]/2.) * sum(armatau2.col(j)))));
   
      // draw tau^2
      for (int k = 0; k < m; k++) {
       if (armarestr(k,j) != 0) {
        armatau2(k,j) = do_rgig1(aShrink(j) - .5, armafacload(k,j)*armafacload(k,j),
	 aShrink(j)*armalambda2(j));
       }
      }
     }
    }
   }
   // should we employ the DL-prior?
   if (dlprior) {
    tmpcounter4samplingtauDL = 0;
    for (int j = 0; j < r; j++) {
     for (int k = 0; k < m; k++) {
      if (armarestr(k,j) != 0) {
       armapsiDL(k,j) = 1. / do_rgig1(-.5, 1, (armafacload(k,j)*armafacload(k,j)) / (tauDL * tauDL * armaphiDL(k,j) * armaphiDL(k,j)));
       tmpcounter4samplingtauDL += fabs(armafacload(k,j))/armaphiDL(k,j);
       armaTDL(k,j) = do_rgig1(aShrink[0] - 1., 2*fabs(armafacload(k,j)), 1);
      }
     }
    }
    //Rprintf("%f\n", tmpcounter4samplingtauDL);
    //if (tmpcounter4samplingtauDL < 1.) Rprintf("THIS: %f\n", armaphiDL[0,0]);
    tauDL = do_rgig1((aShrink[0] - 1.) * unrestrictedelementcount, 2. * tmpcounter4samplingtauDL, 1);
    armaphiDL = armaTDL / accu(armaTDL);
    armatau2 = armapsiDL % armaphiDL % armaphiDL * tauDL * tauDL;
   }
  
   int oldpos = 0;
   for (int j = 0; j < m; j++) {
      
   // TODO: some things outside
   
   
   // transposed design matrix Xt is filled "manually"
   int activecols = 0;
   for (int l = 0; l < r; l++) {
    if (armarestr(j, l) != 0) {
     for (int k = 0; k < T; k++) {
      armaXt(activecols, k) = armaf(l, k) * armahtilde(k, j);
     }
     activecols++;
    }
   }

   armaytilde = armay.row(j).t() % armahtilde.col(j);
  

   // Now draw from the multivariate normal distribution
   // armaSigma is first used as temporary variable:
   armaSigma.submat(0,0,activecols-1,activecols-1) = armaXt.rows(0,activecols-1) * armaXt.rows(0,activecols-1).t();

   // add precisions to diagonal:
   armaSigma.submat(0,0,activecols-1,activecols-1).diag() += 1/arma::nonzeros(armatau2.row(j));

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
   armamean.head(activecols) = armaSigma.submat(0, 0, activecols-1, activecols-1) *
                                armaXt.submat(0, 0, activecols-1, T-1) *
                                armaytilde;
   
   // draw from the r-variate normal distribution
   
   armadraw = rnorm(r);
   
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
  armafacloadt(armafacloadtunrestrictedelements) = armafacloadtmp;
  armafacload = arma::trans(armafacloadt);

 //Rprintf("\n\n");
 //for (int is = 0; is < m; is++) Rprintf("%f %f\n", curfacload(is, 0), curfacload(is, 1));

  if (interweaving == 1 || interweaving == 3 || interweaving == 5 || interweaving == 7) {
  // STEP 2*: "Shallow" Interweaving
  
//   // intermediate step: calculate transformation of curh
//   armahtilde2 = exp(-armah(arma::span::all, arma::span(m, m+r-1)));
   

   for (int j = 0; j < r; j++) {
    
    int userow = j;
    if (interweaving == 3 || interweaving == 7) { // find largest absolute element in column to interweave
     userow = 0;
     for (int k = 1; k < m; k++) if (fabs(armafacload(k, j)) > fabs(armafacload(userow, j))) userow = k;
    }
    if (interweaving == 5) { // find random nonzero element in column to interweave
     for (int k = 1; k < m; k++) {
      userow = floor(R::runif(0, m));
      if (fabs(armafacload(userow, j)) > 0.01) break;
     }
    }

     
    double newdiag2 = do_rgig1((nonzerospercol(j)- T) / 2.,
       sum(square(armaf.row(j) * armafacload(userow,j))/exp(armah.col(m+j)).t()),
       sum(square(nonzeros(armafacload.col(j))) / nonzeros(armatau2.col(j))) / (armafacload(userow,j) * armafacload(userow,j)));
   
    double tmp = sqrt(newdiag2)/armafacload(userow,j);
   
    armafacload.col(j) *= tmp;
    armaf.row(j) *= 1/tmp;
   }
  }

  if (interweaving == 2 || interweaving == 4 || interweaving == 6 || interweaving == 7) { // STEP 2+: "Deep" Interweaving 
   for (int j = 0; j < r; j++) {

    int userow = j;
    if (interweaving == 4 || interweaving == 7) { // find largest absolute element in column to interweave
     userow = 0;
     for (int k = 1; k < m; k++) if (fabs(armafacload(k, j)) > fabs(armafacload(userow, j))) userow = k;
    }
    if (interweaving == 6) { // find random nonzero element in column to interweave
     for (int k = 1; k < m; k++) {
      userow = floor(R::runif(0, m));
      if (fabs(armafacload(userow, j)) > 0.01) break;
     }
      //Rprintf("use: %i ", userow);
    }


    //Rprintf("%i and %i\n", j, userow);

    double phi = curpara(1,m+j); 
    double sigma = curpara(2,m+j);
    double mu_old = log(armafacload(userow,j) * armafacload(userow,j));
    hopen = curh(_, m+j) + mu_old;
    double h0open = curh0(m+j) + mu_old;
    double logacceptrate;
    double mu_prop;

    if (priorh0(m+j) < 0.) {  // old prior for h0 (stationary distribution, depends on phi), as in JCGS submission Feb 2016
     double tmph = hopen(0) - phi*h0open;
     for (int k = 1; k < T; k++) tmph += hopen(k) - phi*hopen(k-1);
   
     double gamma_old = (1 - phi) * mu_old;
     double gamma_prop = as<double>(rnorm(1, tmph/(T+B011inv), sigma/sqrt(T+B011inv)));
     mu_prop = gamma_prop/(1-phi);

     logacceptrate = logdnormquot(mu_prop, mu_old, h0open, sigma/sqrt(1-phi*phi));
     logacceptrate += logspecialquot(gamma_prop, gamma_old, .5, 1/(2.*armatau2(userow,j)), 1-phi);
     logacceptrate += logdnormquot(gamma_old, gamma_prop, 0., sigma*sqrt(1/B011inv));
   
    } else {  // new prior does not depend on phi
     double tmph = hopen(0);
     for (int k = 1; k < (T-1); k++) tmph += hopen(k);
     
     double tmp4prop = T*priorh0(m+j)*(1-phi)*(1-phi) + 1;
     double prop_mean = (priorh0(m+j) * (1-phi) * (hopen(T-1) + (1-phi)*tmph - phi*h0open) + h0open) / tmp4prop;
     double prop_sd = (sqrt(priorh0(m+j)) * sigma) / sqrt(tmp4prop);
    
     mu_prop = as<double>(rnorm(1, prop_mean, prop_sd));
     logacceptrate = .5 * ((mu_prop - mu_old) - (exp(mu_prop) - exp(mu_old)) / armatau2(userow,j));
    }
    
    // NEW, same for both priors:
    arma::vec relevantload = armafacload.col(j);
    arma::vec relevanttau2 = armatau2.col(j);
    
    // use all except interwoven element (restricted loadings are assumed to be zero!)
    double mysum = accu(square(nonzeros(relevantload))/nonzeros(relevanttau2)) -
     (relevantload(userow)*relevantload(userow))/relevanttau2(userow); 
    
    logacceptrate += .5 * ((nonzerospercol(j)-1)*(mu_prop - mu_old) -
     mysum / (armafacload(userow,j)*armafacload(userow,j)) * (exp(mu_prop) - exp(mu_old)));
     

   // Rprintf("ACCEPT? ");

    //ACCEPT/REJECT
    if (log(::unif_rand()) < logacceptrate) {
//    Rprintf("ACC col %i el %02i - ", j+1, userow+1);
     curh(_, m+j) = hopen - mu_prop;
     curh0(m+j) = h0open - mu_prop;

     double tmp = exp(mu_prop/2)/armafacload(userow,j);
     armafacload.col(j) *= tmp;
     armaf.row(j) *= 1/tmp;
//    } else {
//     Rprintf("REJ col %i el %02i - ", j+1, userow+1);
    }
   }
//   Rprintf("\n");
  }
  // STEP 3:
  // update the factors (T independent r-variate regressions with m observations)

  if (samplefac) {
  armadraw2.imbue(::norm_rand);
  for (int j = 0; j < T; j++) {
   
   // transposed design matrix Xt2 (r x m) is filled "manually"
   for (int k = 0; k < m; k++) {
    for (int l = 0; l < r; l++) {
     Xt2(l, k) = armafacload(k, l) * armahtilde(j,k);
    }
   }
   
   armaytilde2 = armay.col(j) % armahtilde.row(j).t();
  
   // Now draw form the multivariate normal distribution

   // armaSigma2 is first used as temporary variable (to hold the precision):
   armaSigma2 = armaXt2 * armaXt2.t();
   
   // add precisions to diagonal:
   armaSigma2.diag() += exp(-armah(j, arma::span(m, mpr-1)));

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
   armamean2 = armaSigma2 * armaXt2 * armaytilde2;

   // draw from the r-variate normal distribution
   try {
    armaf.col(j) = armamean2 + (armaR2inv * armadraw2.subvec(j*r, (j+1)*r - 1));
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
     armafacload.col(j) *= -1;
     armaf.row(j) *= -1;
    }
   }
  }

  // REGRESSION (only single regression)
  if (model_mean) {
   if (r > 0) {
    armay_regression = armay_original - armafacload * armaf;
   }
   for (int j = 0; j < m; j++) {
    const arma::vec expmh2 = arma::exp(-.5*armah.col(j));
    beta_j[0] = armabeta[j];
    stochvol::update_regressors(armay_regression.row(j).t() % expmh2, expmh2, beta_j, prior_specs[j]);
    armabeta[j] = beta_j[0];
   }

  // de-meaned observations
  armay = armay_original;
  armay.each_col() -= armabeta;
  }

  // STORAGE:
  if (effi >= 0) {
   if (effi % thin == (thin - 1)) {
   store(curfacload, facload, curf, f, curh, h, curh0, h0,
     curpara, para, curlambda2, lambda2, curtau2, tau2, curbeta,
     beta, curmixind, mixind, auxstore, thintime, (i-burnin)/thin);
   }
 //Rprintf("\n");
 //for (int is = 0; is < m; is++) Rprintf("%f %f\n", facload(i*m*r+is), facload(i*m*r+m+is));

   if (effi % runningstoreevery == (runningstoreevery - 1)) {
    if (runningstore >= 1) {
     armahrunmean = (armahrunmean * effirunningstore + armah) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armahrunm2 = (armahrunm2 * effirunningstore + pow(armah, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armahrunm3 = (armahrunm3 * effirunningstore + pow(armah, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armahrunm4 = (armahrunm4 * effirunningstore + pow(armah, 4)) / (effirunningstore + 1);
//     armahrunmin = min(armahrunmin, armah);
//     armahrunmax = max(armahrunmax, armah);
    }

    if (runningstore >= 2) {
     armafrunmean = (armafrunmean * effirunningstore + armaf) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armafrunm2 = (armafrunm2 * effirunningstore + pow(armaf, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armafrunm3 = (armafrunm3 * effirunningstore + pow(armaf, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armafrunm4 = (armafrunm4 * effirunningstore + pow(armaf, 4)) / (effirunningstore + 1);
//     armafrunmin = min(armafrunmin, armaf);
//     armafrunmax = max(armafrunmax, armaf);
    }

    if (runningstore >= 3) {
     htranstmp = exp(armah/2);
     armahrunmeantrans = (armahrunmeantrans * effirunningstore + htranstmp) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armahrunm2trans = (armahrunm2trans * effirunningstore + pow(htranstmp, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armahrunm3trans = (armahrunm3trans * effirunningstore + pow(htranstmp, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armahrunm4trans = (armahrunm4trans * effirunningstore + pow(htranstmp, 4)) / (effirunningstore + 1);
    }
   
    if (runningstore >= 4) {
     tmpcov.fill(0.);
     tmpcounter = 0;
     for (int k = 0; k < m; k++) {
      for (int l = k; l < m; l++) {
       for (int inner = 0; inner < r; inner++) {
        tmpcov.col(tmpcounter) += armafacload(l,inner)*armafacload(k,inner)*exp(armah.col(m+inner));
       }
       if (k == l) {
        tmpcov.col(tmpcounter) += exp(armah.col(k));
        tmpvol.col(k) = sqrt(tmpcov.col(tmpcounter));
       }
       tmpcounter++;
      }
     }
     armacovrunmean = (armacovrunmean * effirunningstore + tmpcov) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armacovrunm2 = (armacovrunm2 * effirunningstore + pow(tmpcov, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armacovrunm3 = (armacovrunm3 * effirunningstore + pow(tmpcov, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armacovrunm4 = (armacovrunm4 * effirunningstore + pow(tmpcov, 4)) / (effirunningstore + 1);
     armavolrunmean = (armavolrunmean * effirunningstore + tmpvol) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armavolrunm2 = (armavolrunm2 * effirunningstore + pow(tmpvol, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armavolrunm3 = (armavolrunm3 * effirunningstore + pow(tmpvol, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armavolrunm4 = (armavolrunm4 * effirunningstore + pow(tmpvol, 4)) / (effirunningstore + 1);
    }

    if (runningstore >= 5) {
     tmpsds = sqrt(tmpcov.cols(diagindices));
     tmpcounter = 0;
     for (int k = 0; k < m; k++) {
      for (int l = k + 1; l < m; l++) {
       tmpcor.col(tmpcounter) = tmpcov.col(tmpcounter + k + 1) / (tmpsds.col(k) % tmpsds.col(l));
       tmpcounter++;
      }
     }
     armacorrunmean = (armacorrunmean * effirunningstore + tmpcor) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armacorrunm2 = (armacorrunm2 * effirunningstore + pow(tmpcor, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armacorrunm3 = (armacorrunm3 * effirunningstore + pow(tmpcor, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armacorrunm4 = (armacorrunm4 * effirunningstore + pow(tmpcor, 4)) / (effirunningstore + 1);
    }

    if (runningstore >= 6) {
     curcom.cols(0,m-1) = 1 - htranstmp.cols(0,m-1) % htranstmp.cols(0,m-1) / tmpcov.cols(diagindices);
     curcom.col(m) = sum(curcom.cols(0,m-1), 1)/m;
     armacomrunmean = (armacomrunmean * effirunningstore + curcom) / (effirunningstore + 1);
     if (runningstoremoments >= 2) armacomrunm2 = (armacomrunm2 * effirunningstore + pow(curcom, 2)) / (effirunningstore + 1);
     if (runningstoremoments >= 3) armacomrunm3 = (armacomrunm3 * effirunningstore + pow(curcom, 3)) / (effirunningstore + 1);
     if (runningstoremoments >= 4) armacomrunm4 = (armacomrunm4 * effirunningstore + pow(curcom, 4)) / (effirunningstore + 1);
    }
    
    effirunningstore += 1;
   }
  }
 }  // END main MCMC loop

 if (verbose) {
   Rprintf("\r********* Iteration %*i of %*i (%3.0f%%) *********",
    space4print, N, space4print, N, 100.);
 }
 
 List retval = List::create(
   Named("facload") = facload,
   Named("fac") = f,
   Named("logvar") = h,
   Named("logvar0") = h0,
   Named("para") = para,
   Named("beta") = beta,
   Named("mixind") = mixind,
   Named("lambda2") = lambda2,
   Named("tau2") = tau2,
   Named("latestauxiliary") = List::create(
     Named("lambda2") = curlambda2,
     Named("facloadvar") = curtau2),
   Named("y") = y,
   Named("runningstore") = List::create(
     Named("logvar") = List::create(
       Named("mean") = hrunmean,
       Named("m2") = hrunm2,
       Named("m3") = hrunm3,
       Named("m4") = hrunm4),
     Named("fac") = List::create(
       Named("mean") = frunmean,
       Named("m2") = frunm2,
       Named("m3") = frunm3,
       Named("m4") = frunm4),
     Named("sd") = List::create(
       Named("mean") = hrunmeantrans,
       Named("m2") = hrunm2trans,
       Named("m3") = hrunm3trans,
       Named("m4") = hrunm4trans),
     Named("cov") = List::create(
       Named("mean") = covrunmean,
       Named("m2") = covrunm2,
       Named("m3") = covrunm3,
       Named("m4") = covrunm4),
     Named("vol") = List::create(
       Named("mean") = volrunmean,
       Named("m2") = volrunm2,
       Named("m3") = volrunm3,
       Named("m4") = volrunm4),
     Named("cor") = List::create(
       Named("mean") = corrunmean,
       Named("m2") = corrunm2,
       Named("m3") = corrunm3,
       Named("m4") = corrunm4),
     Named("com") = List::create(
       Named("mean") = comrunmean,
       Named("m2") = comrunm2,
       Named("m3") = comrunm3,
       Named("m4") = comrunm4)));
 //PutRNGstate();
 return retval;
}
