#include "sampler.h"
#include "predict.h"
#include "dmvnorm.h"

using namespace Rcpp;

static const R_CallMethodDef CallEntries[] = {
    {"sampler", (DL_FUNC) &sampler, 32},
    {"predict", (DL_FUNC) &predict, 3},
    {"dmvnorm", (DL_FUNC) &dmvnorm, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_factorstochvol(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
