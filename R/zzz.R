#  #####################################################################################
#  R package factorstochvol by
#     Gregor Kastner Copyright (C) 2016-2020
#     Darjus Hosszejni Copyright (C) 2019-2020
#  
#  This file is part of the R package factorstochvol: Bayesian Estimation
#  of (Sparse) Latent Factor Stochastic Volatility Models
#  
#  The R package factorstochvol is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 2 or any
#  later version of the License.
#  
#  The R package factorstochvol is distributed in the hope that it will
#  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with the R package factorstochvol. If that is not the case,
#  please refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################

.onAttach <- function(lib, pkg) {
	if(interactive() || getOption("verbose")){
		packageStartupMessage(sprintf("Package %s %s attached. To cite, see citation(\"%s\").", pkg, utils::packageDescription(pkg)$Version, pkg))
	}
}

.onUnload <- function (libpath) {
  library.dynam.unload("factorstochvol", libpath)
}
