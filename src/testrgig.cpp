#include "sampler.h"

using namespace Rcpp;

RcppExport SEXP testrgig() {
 double tmp = rgig1(1.,2.,3.);
 return List::create(Named("value") = tmp);
}
