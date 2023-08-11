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

// Main sampler (as called from R):
RcppExport SEXP sampler(const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP);

#endif