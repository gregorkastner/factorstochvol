#ifndef RCPPEXPORTS_H
#define RCPPEXPORTS_H

#include <RcppArmadillo.h>

RcppExport SEXP _factorstochvol_dmvnorm(SEXP x_SEXP, SEXP means_SEXP, SEXP vars_SEXP, SEXP logaSEXP);

RcppExport SEXP _factorstochvol_predict(SEXP objSEXP, SEXP storeSEXP, SEXP eachSEXP);
  
RcppExport SEXP _factorstochvol_sampler(SEXP ySEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP startval_inSEXP, SEXP bmuSEXP, SEXP BmuSEXP, SEXP priorphiSEXP, SEXP BsigmaSEXP, SEXP priorbetaSEXP, SEXP model_meanSEXP, SEXP shrinkagepriorsSEXP, SEXP thinSEXP, SEXP auxstoreSEXP, SEXP thintimeSEXP, SEXP quietSEXP, SEXP parameterizationSEXP, SEXP MHstepsSEXP, SEXP B011_inSEXP, SEXP B022_inSEXP, SEXP MHcontrolSEXP, SEXP GammapriorSEXP, SEXP offsetSEXP, SEXP truncnormalSEXP, SEXP restrSEXP, SEXP interweavingSEXP, SEXP signswitchSEXP, SEXP runningstoreSEXP, SEXP runningstoreeverySEXP, SEXP runningstoremomentsSEXP, SEXP pflSEXP, SEXP svSEXP, SEXP priorhomoskedasticSEXP, SEXP priorh0SEXP, SEXP samplefacSEXP, SEXP facloadtolSEXP);

#endif