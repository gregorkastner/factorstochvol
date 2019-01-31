#ifndef _PREDICT_H
#define _PREDICT_H

//#define ARMA_NO_DEBUG // disables bounds checks
#include <RcppArmadillo.h>

// Main predict function (as called from R):
RcppExport SEXP predict(const SEXP, const SEXP, const SEXP);

#endif
