#  #####################################################################################
#  R package factorstochvol by
#     Gregor Kastner Copyright (C) 2016-2021
#     Darjus Hosszejni Copyright (C) 2019-2021
#     Luis Gruber Copyright (C) 2021
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

#' Plot communalities over time.
#'
#' \code{comtimeplot} plots the communalities over time, i.e. the
#' series-specific percentage of variance explained through the common factors.
#'
#' This function displays the joint (average) communalities over time and all
#' series-specific communalities. If communalities haven't been stored during
#' sampling, \code{comtimeplot} produces an error.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj Object of class \code{'fsvsim'} (or \code{NULL}), usually
#' resulting from a call to \code{\link{fsvsim}}. Defaults to \code{NULL}.
#' @param show Indicator whether to show joint (\code{'joint'}), series-specific
#' (\code{'series'}), or both (\code{'both'}) communalities.
#' @param maxrows Single positive integer denoting the maximum number of series
#' in each plot. Defaults to 5.
#' @param ylim Vector of length two denoting the range of the horizontal axis.
#' Defaults to 1.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @import methods
#'
#' @export
comtimeplot <- function(x, fsvsimobj = NULL, show = "series",
			maxrows = 5, ylim = c(0,100)) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!exists("runningstore", x) || !exists("com", x$runningstore))
  stop("Communalities haven't been stored during sampling.")
 if (length(show) != 1) stop("Illegal value of 'show'.")

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]
 snames <- colnames(x$y)
 dates <- rownames(x$y)
 if (is.null(dates)) dates <- 1:n

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
 }

 if (!is.null(fsvsimobj)) {
  communalities <- 1 - exp(t(fsvsimobj$idivol)) / apply(covmat(fsvsimobj),3,diag)
 }

 for (i in c(m+1,1:m)) {
  if (i == m+1 && show == "series") next
  if (i < m+1 && show == "joint") next
  if (i == 1) oldpar <- par(mfrow = c(min(maxrows, m), 1),
			    mgp = c(2, .5, 0),
			    mar = c(1.5, 1.5, 1.5, 0.5))
  thismean <- x$runningstore$com[,i,"mean"]
  thissd <- x$runningstore$com[,i,"sd"]
  ts.plot(100*cbind(thismean - 2*thissd, thismean, thismean + 2*thissd),
	  lwd = c(1, 2, 1), main = "", xlab = "", ylab = "",
	  gpars = list(xaxt = 'n'), ylim = ylim)
  ats <- seq(1, n, len = 11)
  axis(1, labels = dates[ats], at = ats)

  if (i == m+1) {
   title(paste0("Joint communalities (mean +/- 2sd)"))
   if (!is.null(fsvsimobj)) lines(colMeans(communalities), col = 3)
  } else {
   title(paste0("Communalities of series ", i, " (", snames[i], ", mean +/- 2sd)"), font.main = 1)
   if (!is.null(fsvsimobj)) lines(communalities[i,], col = 3)
  }
 }
 if (show != "joint") par(oldpar)
 invisible(x)
}


#' Plot series-specific volatilities over time.
#'
#' \code{voltimeplot} plots the marginal volatilities over time, i.e. the
#' series-specific conditional standard deviations. If these haven't been
#' stored during sampling (because \code{runningstore} has been set too low),
#' \code{voltimeplot} throws a warning.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param these Index vector containing the time points to plot. Defaults
#' to \code{seq_len(nrow(x$y))}, i.e., all timepoints.
#' @param legend Where to position the \code{\link{legend}}.
#' If set to NULL, labels will be put directly next to the series.
#' Defaults to "topright".
#' @param ... Additional parameters will be passed on to \code{\link{ts.plot}}.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
voltimeplot <- function(x, these = seq_len(nrow(x$y)), legend = "topright", ...) {
if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!exists("runningstore", x) || !exists("vol", x$runningstore)) {
  warning("Implied volatilities haven't been stored during sampling.")
  return(invisible(x))
 }
 n <- nrow(x$y)
 m <- ncol(x$y)
 r <- dim(x$logvar)[2] - m
 draws <- dim(x$para)[3]
 snames <- colnames(x$y)
 dates <- rownames(x$y)

 if (!is.numeric(these) || min(these) < 1 || max(these) > n) stop("Illegal argument value 'these'.")

 dat <- matrix(x$runningstore$vol[these,,"mean",drop = FALSE], nrow = length(these))
# plotorder <- order(colMeans(dat))
 plotorder <- seq_len(ncol(dat))
 if (length(palette()) != m) colas <- rainbow(m) else colas <- seq_len(m)

 ts.plot(dat[,plotorder,drop = FALSE], gpars = list(col = colas, xaxt = 'n', xlab = '', ...))
 ats <- round(seq(1, length(these), length.out = min(length(these), 10)))
 axis(1, at = ats, labels = dates[these][ats])

 mynames <- snames[plotorder]
 if (!is.null(legend)) {
  legend(legend, legend = mynames, col = colas, ncol = 2, lty = 1, lwd = 2)
 } else {
  text(-.018*length(these), dat[1,], mynames, col = colas)
  text(1.018*length(these), dat[nrow(dat),], mynames, col = colas)
 }
 invisible(x)
}

#' Plot correlation matrices for certain points in time
#'
#' \code{corimageplot} plots the model-implied correlation matrices
#' for one or several points in time.
#'
#' @note If correlations haven't been stored during sampling,
#' \code{corimageplot} produces an error.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param these Index vector containing the time points to plot. Defaults
#' to \code{seq_len(nrow(x$y))}.
#' @param order String, where \code{'none'} and \code{'original'}
#' indicate not to mess
#' with the series ordering. Other keywords
#' (e.g. \code{'hclust'}) will be forwarded to
#' \code{\link[corrplot]{corrMatOrder}}.
#' @param these4order Index vector containing the time points used for
#' ordering. Probably, the default (\code{these}) is what you want.
#' @param plotdatedist Numerical value indicating where the dates should
#' be plotted.
#' @param plotCI String. If not equal to 'n', posterior credible regions are
#' added (posterior mean +/- 2 posterior sd). Ignored if \code{plottype} is
#' "imageplot".
#' @param date.cex Size multiplier for the dates.
#' @param col Color palette or NULL (the default).
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param plottype Indicates which type of plot should be drawn. Can be
#' "corrplot" for \code{\link[corrplot]{corrplot}} (recommended for up to around
#' 20 series), or "imageplot" for a simpler \code{\link[graphics]{image}} plot.
#' @param ... Additional parameters will be passed on to
#' \code{\link[corrplot]{corrplot}}. Ignored if \code{plottype} is
#' "imageplot".
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
corimageplot <- function(x, these = seq_len(nrow(x$y)), order = "original",
			 these4order = these,
			 plotdatedist = 0, plotCI = 'n', date.cex = 1.5, col = NULL,
			 fsvsimobj = NULL, plottype = "corrplot", ...) {
 type <- "cor"
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }
 if (!exists("runningstore", x) || !exists(type, x$runningstore))
  stop("Correlations haven't been stored during sampling.")
 n <- nrow(x$y)
 m <- ncol(x$y)
 r <- dim(x$logvar)[2] - m
 draws <- dim(x$para)[3]
 snames <- colnames(x$y)
 if (is.null(snames)) snames <- 1:m
 dates <- rownames(x$y)
 if (is.null(dates)) dates <- 1:n

 if (!is.numeric(these) || min(these) < 1 || max(these) > n) stop("Illegal argument value 'these'.")

 if (order != 'none' && order != 'original') {
  orderthis <- matrix(NA_real_, nrow = length(these4order), ncol = m)
  for (i in seq(along = these4order)) {
   this <- runningcormat(x, these4order[i], type = type, statistic = "mean")
   orderthis[i,] <- corrplot::corrMatOrder(this, order = order)
  }
 } else orderthis <- matrix(1:m, nrow = length(these), ncol = m, byrow = TRUE)
 tmp <- apply(orderthis, 1, paste, collapse = ' ')
 orderthis <- names(which.max(table(tmp)))
 orderthis <- as.integer(unlist(strsplit(orderthis, " ")))

 if (is.null(col)) {
  colpal <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
    "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  col <- rev(colpal(200))
 }

 for (i in these) {
  toplot <- runningcormat(x, i, type = type, statistic = "mean")
  rownames(toplot) <- colnames(toplot) <- snames
  if (!(plotCI == 'n')) {
   toplot2 <- runningcormat(x, i, type = type, statistic = "sd")
   lower <- pmax(toplot - 2 * toplot2, -1)
   upper <- pmin(toplot + 2 * toplot2, 1)
   rownames(lower) <- colnames(lower) <- snames
   rownames(upper) <- colnames(upper) <- snames
  } else lower <- upper <- NULL

  if (plottype == "corrplot") {
   corrplot::corrplot(toplot[orderthis, orderthis], plotCI = plotCI,
		      lowCI.mat = lower[orderthis, orderthis],
		      uppCI.mat = upper[orderthis, orderthis], diag = FALSE,
		      col = col, ...)
   text((m+1)/2, m+plotdatedist, dates[i], cex = date.cex)
  } else if (plottype == "imageplot") {
   image(toplot[orderthis, orderthis])
  }

  if (!is.null(fsvsimobj)) {
   cortrue <- cov2cor(covmat(fsvsimobj, i)[,,1])
   diag(cortrue) <- NA
   oldpar <- par(xpd = TRUE)
   symbols(rep(1:m, each = m), rep(m:1, m),
	   circles = .9*abs(as.numeric(cortrue[orderthis, orderthis]))^0.5/2,
	   fg = "green", inches = FALSE, add = TRUE)
   par(oldpar)
  }
 }

 invisible(x)
}


#' Plot correlations over time.
#'
#' \code{cortimeplot} draws correlations over time.
#'
#' This function displays one component series' time-varying correlations with
#' the other components series. Throws an error if correlations haven't been
#' stored during sampling.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param series Single number, coercible to integer. Indicates the series
#' relative to which correlations are drawn.
#' @param these Index vector containing the time points to plot. Defaults
#' to \code{seq_len(nrow(x$y))}.
#' @param type What to plot, usually "cor" or "cov".
#' @param statistic Which posterior summary should be plotted, usually "mean".
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
cortimeplot <- function(x, series, these = seq_len(nrow(x$y)),
			type = "cor", statistic = "mean") {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!exists("runningstore", x) || !exists(type, x$runningstore))
  stop("What you are requesting (argument 'type') hasn't been stored during sampling.")
 if (missing(series) || !is.numeric(series) || length(series) != 1 || series < 1) {
  stop("Argument 'series' must be a single number, coercible to integer.")
 } else series <- as.integer(series)
 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$logvar)[2] - m
 draws <- dim(x$para)[3]
 snames <- colnames(x$y)
 dates <- rownames(x$y)
 if (is.null(dates)) dates <- 1:n
 oldpar <- par(mgp = c(1.8, .6, 0), mar = c(1.7, 1.7, 1.7, .3))

 if (!is.numeric(these) || length(these) < 1 || any(these > n) || any(these < 1)) {
  stop("Argument 'these' must be numeric with elements between 1 and nrow(x$y).")
 }

 tryCatch(tmp <- x$runningstore[[type]], error = function(e)
	   stop("Argument 'type' must be either 'cov' or 'cor'."))

 tryCatch(meancors <- tmp[these,,statistic], error = function(e)
	   stop(paste0("Argument 'statistic' must be one of: ",
            paste(dimnames(tmp)[[3]], collapse = ', '), ".")))

 curind1 <- grep(paste0('^', series, '_'), colnames(meancors))
 curind2 <- grep(paste0('_', series, '$'), colnames(meancors))

 if (type == "cor") curind <- c(curind1, curind2)
 if (type == "cov") curind <- curind2

 colororder <- order(colMeans(meancors[,curind]))
 curcors <- meancors[,curind[colororder]]

 if (type == "cor") cornames <- snames[-series][colororder]
 if (type == "cov") cornames <- snames[colororder]


 if (length(palette()) != m) colas <- rainbow(m) else colas <- seq_len(m)

 ts.plot(curcors, col = colas, gpars = list(xaxt = 'n', xlab = '', ylab = ''))
 title(paste0('Posterior ', statistic, ' of pairwise ', type, 's with ', snames[series]))
 abline(h = 0, lty = 3)
 myseq <- round(seq(1, length(these), length.out = 10))
 axis(1, labels = dates[these][myseq], at = myseq)
 text(-.018*length(these), curcors[1,], cornames, col = colas)
 text(1.018*length(these), curcors[nrow(curcors),], cornames, col = colas)

 par(oldpar)

 invisible(x)
}


#' @rdname cortimeplot
#' @export
covtimeplot <- function(x, series, these = seq_len(nrow(x$y)),
			type = "cov", statistic = "mean") {
 cortimeplot(x = x, series = series, these = these, type = type,
	     statistic = statistic)
}


#' Displays bivariate marginal posterior distributions of factor loadings.
#'
#' \code{facloadpairplot} illustrates the bivariate marginals of the
#' factor loadings distribution. For a monochrome variant, see
#' \code{\link{facloadcredplot}}.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param maxpoints The maximum amount of posterior draws to plot. If
#' the number of draws stored in \code{x} exceeds this number, draws are
#' thinned accordingly.
#' @param alpha Level of transparency.
#' @param cex Controls the size of the dots.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
facloadpairplot <- function(x, maxpoints = 500, alpha = 20/maxpoints, cex = 3) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (any(dim(x$facload) < 2)) stop("Currently implemented for two or more factors.")

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]
 plotthese <- sample.int(draws, min(draws, maxpoints))

 means <- apply(x$facload, 1:2, mean)

 whiches <- matrix((1:(2*ceiling(r/2))-1) %% r + 1, nrow = 2)

 alpha <- min(alpha, 1)
 if (length(palette()) != m) oldpal <- palette(rainbow(m))

 colas <- apply(sapply(palette(), col2rgb)/255, 2,
                function (cc, alpha) rgb(cc[1], cc[2], cc[3], alpha = alpha), alpha)

 for (i in seq.int(ncol(whiches))) {
  tmp <- aperm(x$facload[,whiches[,i],plotthese], c(1,3,2))
  dim(tmp) <- c(m*length(plotthese), 2)
  myxlims <- quantile(x$facload[,whiches[1,i],plotthese], c(.1/m, 1-.1/m))
  myylims <- quantile(x$facload[,whiches[2,i],plotthese], c(.1/m, 1-.1/m))
  plot(tmp, pch = 16, cex = cex, col = colas, xlim = myxlims, ylim = myylims,
       xlab = paste("Loadings on Factor", whiches[1,i]),
       ylab = paste("Loadings on Factor", whiches[2,i]))
  text(means[,whiches[2*i-1]], means[,whiches[2*i]], colnames(x$y))
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
 }

 if (exists("oldpal")) palette(oldpal)

 invisible(x)
}

#' Displays bivariate marginal posterior distribution of factor loadings.
#'
#' \code{facloadcredplot} illustrates the bivariate marginals of the
#' factor loadings distribution. It is a monochrome variant of
#' \code{\link{facloadpairplot}}.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param quants Posterior quantiles to be plotted.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
facloadcredplot <- function(x, quants = c(.01, .99)) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (any(dim(x$facload) < 2)) stop("Currently implemented for two or more factors.")

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]

 whiches <- matrix((1:(2*ceiling(r/2))-1) %% r + 1, nrow = 2)

 lower <- apply(x$facload, 1:2, quantile, quants[1])
 upper <- apply(x$facload, 1:2, quantile, quants[2])
 med <- apply(x$facload, 1:2, median)

 for (i in seq.int(ncol(whiches))) {
  myxlim <- range(lower[,whiches[1,i]], upper[,whiches[1,i]])
  myylim <- range(lower[,whiches[2,i]], upper[,whiches[2,i]])
  xbreak <- .02*diff(myxlim)
  ybreak <- .02*diff(myylim)
  plot(med[,whiches[,i]], pch = "",
       xlab = paste("Loadings on Factor", whiches[1,i]),
       ylab = paste("Loadings on Factor", whiches[2,i]),
       xlim = myxlim, ylim = myylim)
  for (j in 1:m) {
   lines(c(lower[j,whiches[1,i]], med[j,whiches[1,i]] - xbreak), rep(med[j,whiches[2,i]], 2))
   lines(c(med[j,whiches[1,i]] + xbreak, upper[j,whiches[1,i]]), rep(med[j,whiches[2,i]], 2))

   lines(rep(med[j,whiches[1,i]], 2), c(lower[j,whiches[2,i]], med[j,whiches[2,i]] - ybreak))
   lines(rep(med[j,whiches[1,i]], 2), c(med[j,whiches[2,i]] + ybreak, upper[j,whiches[2,i]]))
  }
  text(med[,whiches[1,i]], med[, whiches[2,i]], colnames(x$y))
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
 }
 invisible(x)
}

#' Displays point estimates of the factor loadings posterior.
#'
#' \code{facloadpointplot} illustrates point estimates (mean, median, ...)
#' of the estimated factor loadings matrix.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param statistic Character string indicating which posterior statistic
#' should be displayed.
#' @param cex Controls the size of the dots.
#' @param alpha Controls the level of transparency.
#' @param allpairs Logical value; if set to TRUE, all possible
#' pairwise combinations will be plotted.
#' @param col Vector of length \code{m} (number of component series),
#' containing \code{\link[grDevices]{rgb}}-type color codes used for
#' plotting. Will be recycled if necessary.
#'
#' @return Returns \code{x} invisibly, throws a warning if there aren't any
#' factors to plot.
#'
#' @family plotting
#'
#' @export
facloadpointplot <- function(x, fsvsimobj = NULL, statistic = "median",
			     cex = 6.5, alpha = 0.2, allpairs = FALSE, col = NULL) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (any(dim(x$facload) == 0)) {
  warning("There aren't any factor loadings to plot.")
  invisible(x)
 }

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]

 if (is.null(allpairs)) {
  if (r > 10) allpairs <- FALSE else allpairs <- TRUE
 }

 snames <- colnames(x$y)
 if (is.null(snames)) snames <- 1:m

 if (r == 1) {
  facloads <- apply(x$facload, 1, eval(statistic))
 } else {
  facloads <- apply(x$facload, 1:2, eval(statistic))
 }

 myplotit <- function(whiches, facloads, snames, cex = 6.5, transparency = 0.2, col = NULL, fsvsimobj = NULL) {
  basecolor <- rgb(0, 0, 1, transparency)
  if (ncol(facloads) > 1 & facloads[1,2] == 0) leadcolor <- c(1,2) else leadcolor <- basecolor
  simcolor <- rgb(0, 1, 0, transparency)
  m <- nrow(facloads)
  basecolas <- rep(basecolor, m)
  if (!is.null(col)) basecolor <- basecolas <- rep(col, length.out = m)
  basepchs <- rep(21, m)
  if (!is.null(fsvsimobj)) {
   rtrue <- ncol(fsvsimobj$facload)
   if (whiches[1] <= rtrue) {
    xlim <- range(facloads[, whiches[1]], fsvsimobj$facload[, whiches[1]], 0)
   } else {
    xlim <- range(facloads[, whiches[1]], 0)
   }
   if (whiches[2] <= rtrue) {
    ylim <- range(facloads[, whiches[2]], fsvsimobj$facload[, whiches[2]], 0)
   } else {
    ylim <- range(facloads[, whiches[2]], 0)
   }
  } else {
   xlim <- range(facloads[,whiches[1]], 0)
   ylim <- range(facloads[,whiches[2]], 0)
  }
  colas <- basecolas
  colas[whiches] <- leadcolor
  pchs <- basepchs
  pchs[x$identifier[whiches,1]] <- c(24,25)
  plot(facloads[,whiches], pch = pchs, cex = cex, bg = basecolor, col = colas,
       main = paste("Posterior", statistic, "of factor loadings"),
       xlab = paste("Factor", whiches[1]), ylab = paste("Factor", whiches[2]),
       xlim = xlim, ylim = ylim)
  text(facloads[,whiches], snames)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  if (!is.null(fsvsimobj)) {
   tmpfacload <- matrix(0, nrow = nrow(facloads), ncol = ncol(facloads))
   tmpfacload[,1:(min(r, rtrue))] <- fsvsimobj$facload[,1:(min(r, rtrue))]
   points(tmpfacload[,whiches], pch = 21, bg = simcolor, col = simcolor,
	  cex = cex)
   simsnames <- colnames(fsvsimobj$y)
   if (is.null(snames)) simsnames <- 1:m
   text(tmpfacload[,whiches], simsnames, col = 3)
  }
 }

# oldpar <- par(mgp = c(1.7, .5, 0), mar = c(2.7, 2.7, 2, 0.5))

 if (r == 1)  {
  myplot <- barplot(facloads, main = paste("Posterior", statistic, "of factor loadings"),
		    ylab = "Loadings", xlab = "Component series", names = 1:m)
  if (!is.null(fsvsimobj)) {
   points(myplot, fsvsimobj$facload[,1], col = 3, pch = 16, cex = 2)
  }
# } else if (r == 3) {  # TODO
#  mainlab <- "Median factor loadings"
#  cols <- colorspace::diverge_hcl(100)
#  thirddim <- facloads[,3]/max(abs(facloads[,3])) # rescaling to -1:1
#  thirddim <- round(thirddim*(99/2)+50.5)
#  par(mar=c(3.5,3.5,1.6,5.5))
#  plot(facloads[,1:2], pch=21, bg=cols[thirddim], col=cols[thirddim], cex=5.6, main='', xlab='', ylab='')
#  mtext(mainlab, cex=1.35, line=.5)
#  mtext("Factor 1", side=1, cex=1, line=2.5)
#  mtext("Factor 2", side=2, cex=1, line=2.5)
#  text(facloads[,1:2], abbrev)
#  xdist <- diff(range(facloads[,1]))
#  xmax <- max(facloads[,1])
#  ymax <- max(facloads[,2])
#  ymin <- min(facloads[,2])
#  zdist <- 2*ceiling(5*max(abs(facloads[,3])))/5
#  ypos <- seq(ymin, ymax, len=lala <- 5*zdist+2)
#  rect(xmax + .065*xdist, ypos[-lala], xmax + .11*xdist, ypos[-1], col=cols[seq(1,100,len=lala-1)], xpd=NA, border="light grey")
#  text(xmax + .16*xdist, (ypos[-1]+ypos[-lala])/2, round(seq(-zdist/2,zdist/2,len=5*zdist+1), 1),xpd=NA, pos=2)
#  mtext("Factor 3", side=4)
 } else {
  if (allpairs) {
   apply(combn(1:r, 2), 2, myplotit,
	 facloads = facloads, fsvsimobj = fsvsimobj, snames = snames,
	 cex = cex, transparency = alpha, col = col)
  } else {
   apply(matrix((1:(2*ceiling(r/2))-1) %% r + 1, nrow = 2), 2, myplotit,
	 facloads = facloads, fsvsimobj = fsvsimobj, snames = snames,
	 cex = cex, transparency = alpha, col = col)
  }
 }
# par(oldpar)
 invisible(x)
}

#' Plot log-variances over time.
#'
#' \code{logvartimeplot} plots the idiosyncratic and factor log-variances over time.
#'
#' This function displays the posterior distribution (\code{mean +/- 2sd})
#' of log-variances of both
#' the factors and the idiosyncratic series.
#' If these haven't been stored during
#' sampling, \code{logvartimeplot} produces an error.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param show If set to "fac", only factor log-volatilities will be displayed.
#' If set to "idi", only idiosyncratic log-volatilities will be displayed.
#' If set to "both", factor log-volatilities will be drawn first, followed
#' by the idiosyncratic log-volatilities.
#' @param maxrows Indicates the maximum number of rows to be drawn per page.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
logvartimeplot <- function(x, fsvsimobj = NULL, show = "both", maxrows = 5) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!exists("runningstore", x)) stop("Cannot plot time-varying volatilities because
				      'runningstore' was set too low during sampling.")
 if (!is.character(show) || length(show) != 1 || !show %in% c("both", "fac", "idi"))
  stop("Argument 'show' must be one of: 'both', 'fac', 'idi'.")

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]
 snames <- colnames(x$y)
 dates <- rownames(x$y)
 if (is.null(dates)) dates <- 1:n

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }

 oldpar <- par(mfrow = c(min(maxrows, max(r, 1)), 1), mgp = c(2, .5, 0), mar = c(1.5, 1.5, 1.5, 0.5))

 if (show == "both") thesei <- c(m + seq_len(r), seq_len(m))
 else if (show == "fac") thesei <- m + seq_len(r)
 else if (show == "idi") thesei <- seq_len(m)

 for (i in thesei) {
  if (i == 1) par(mfrow = c(min(maxrows, m), 1))
  thismean <- x$runningstore$logvar[,i,"mean"]
  thissd <- x$runningstore$logvar[,i,"sd"]
  ts.plot(cbind(thismean - 2*thissd, thismean, thismean + 2*thissd),
	  lwd = c(1, 2, 1), main = "", xlab = "", ylab = "",
	  gpars = list(xaxt = 'n'))
  ats <- seq(1, n, len = 11)
  axis(1, labels = dates[ats], at = ats)

  if (i <= m) {
   title(paste0("Idiosyncratic log-variance of series ", i, " (", snames[i], ", mean +/- 2sd)"), font.main = 1)
  } else {
   title(paste("Log-variance of factor", i-m, "(mean +/- 2sd)"), font.main = 1)
  }
  if (!is.null(fsvsimobj)) {
   if (i <= m) {
    lines(fsvsimobj$idivol[,i], col = 3)
   } else if (i <= m + rtrue) {
    lines(fsvsimobj$facvol[,i-m], col = 3)
   } else {
    abline(h = 0, col = 2, lty = 3)
   }
  }
 }
 par(oldpar)
 invisible(x)
}

#' Trace plots of parameter draws.
#'
#' \code{paratraceplot} draws trace plots of all parameters (\code{mu, phi,
#' sigma}). Can be an important tool to check MCMC convergence if inference
#' about (certain) parameters is sought.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param thinning Plot every \code{thinning}th draw.
#' @param maxrows Indicates the maximum number of rows to be drawn per page.
#' @param ... Ignored.
#' @return Returns \code{x} invisibly.
#' @family plotting
#' @rdname paratraceplot
#' @name paratraceplot
#' @export
paratraceplot.fsvdraws <- function(x, fsvsimobj = NULL, thinning = NULL, maxrows = 3, ...) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]
 if (is.null(thinning)) thinning <- ceiling(draws/10000)
 plotwhich <- seq(1, draws, by = thinning)
 snames <- colnames(x$y)

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }

 oldpar <- par(mgp = c(2, .5, 0), mar = c(0, 7.4, 0, 2))
 effrows <- min(r, maxrows)
 layout(matrix(1:(3*effrows+2), nrow = 3*effrows+2), heights = c(.075, rep(c(.40, .40, .10)/effrows, effrows), .025))

 for (j in seq_len(r)) {
  if (j %% effrows == 1 | effrows == 1) {
   plot.new()
   text(.5, .7, paste("Parameter draws of factor log-variances with plotthin =", thinning), cex = 1.5, xpd = TRUE)
  }

  plot(x$para[2,m+j,plotwhich], type = "l", main = "", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
  if (!is.null(fsvsimobj)) {
   if(j <= rtrue) {
    abline(h = fsvsimobj$facpara[j,1], col = 3, lty = 2)
   } else {
    legend("topleft", "true value DNE")
   }
  }

  if (j %% effrows == 1) axis(3)
  axis(4)
  axis(2, labels = FALSE)
  mtext(bquote(phi[.(j+m)]), 2, las = 2, line = 2)
  mtext(paste0("Fac", j), 2, las = 2, line = 3.7, padj = (14/effrows), xpd = TRUE)

  plot(x$para[3,m+j,plotwhich], type = "l", main = "", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
  if (!is.null(fsvsimobj)) {
   if (j <= rtrue) {
    abline(h = fsvsimobj$facpara[j,2], col = 3, lty = 2)
   } else {
    legend("topleft", "true value DNE")
   }
  }

  axis(2)
  axis(4, labels = FALSE)
  mtext(bquote(sigma[.(j+m)]), 2, las = 2, line = 2)
  if (j %% effrows == 0) {
   axis(1)
   plot.new()
  }
  plot.new()
 }

 effrows <- min(m, maxrows)
 layout(matrix(1:(4*effrows+2), nrow = 4*effrows+2), heights = c(.075, rep(c(.27, .27, .27, .09)/effrows, effrows), .025))

 for (j in 1:m) {
  if (j %% effrows == 1 | effrows == 1) {
   plot.new()
   text(.5, .7, paste("Parameter draws of idiosyncratic log-variances with plotthin =", thinning), cex = 1.5, xpd = TRUE)
  }

  plot(x$para[1,j,plotwhich], type = "l", main = "", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
  if (!is.null(fsvsimobj)) {
   abline(h = fsvsimobj$idipara[j,1], col = 3, lty = 2)
  }

  if (j %% effrows == 1) axis(3)
  axis(2)
  mtext(bquote(mu[.(j)]), 2, las = 2, line = 1.7)

  plot(x$para[2,j,plotwhich], type = "l", main = "", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
  if (!is.null(fsvsimobj)) {
   abline(h = fsvsimobj$idipara[j,2], col = 3, lty = 2)
  }
  axis(4)
  mtext(bquote(phi[.(j)]), 2, las = 2, line = 1.7)
  mtext(snames[j], 2, las = 2, line = 3.7)

  plot(x$para[3,j,plotwhich], type = "l", main = "", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
  if (!is.null(fsvsimobj)) {
   abline(h = fsvsimobj$idipara[j,3], col = 3, lty = 2)
  }

  axis(2)
  mtext(bquote(sigma[.(j)]), 2, las = 2, line = 1.7)

  if (j %% effrows == 0) {
   axis(1)
   plot.new()
  }

  plot.new()
 }

 par(oldpar)
 invisible(x)
}

#' @export paratraceplot
NULL

#' Trace plots of factor loadings draws
#'
#' \code{facloadtraceplot} draws trace plots of the factor loadings. Can be
#' an important tool to check MCMC convergence if inference about (certain)
#' factor loadings sought.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param thinning Plot every \code{thinning}th draw.
#' @param maxrows Indicates the maximum number of rows to be drawn per page.
#' @param ylim Vector of length two containing lower and upper bounds of the
#' vertical axis. If \code{NULL}, these are automatically determined.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
facloadtraceplot <- function(x, fsvsimobj = NULL, thinning = NULL, maxrows = 10, ylim = NULL) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (any(dim(x$facload) == 0)) {
  warning("There aren't any factor loadings to plot.")
  return(invisible(x))
 }
 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]
 if (is.null(thinning)) thinning <- ceiling(draws/10000)
 plotwhich <- seq(1, draws, by = thinning)
 snames <- colnames(x$y)

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }

 oldpar <- par(mgp = c(2, .5, 0), mar = c(0.1, 6, 0.1, 2))
 effrows <- min(r, maxrows)
 layout(matrix(1:(effrows+2), nrow = effrows+2), heights = c(.075, rep(.9/effrows, effrows), .025))

 for (i in 1:m) {
  for (j in 1:r) {
   if (j %% effrows == 1 | effrows == 1) {
    plot.new()
    text(.5, .5, paste0("Series ", i, " (", snames[i], ") loadings with plotthin = ", thinning), cex = 1.5, xpd = TRUE)
   }
   plot(x$facload[i,j,plotwhich], type = "l", main = "", xlab = "", ylab = "", ylim = ylim, yaxt = 'n', xaxt = 'n')
   if (j %% effrows == 1) axis(3)
   axis(here <- 2 * (j %% 2 + 1))
   axis(6 - here, labels = FALSE)
   mtext(paste("Fac", j), 2, las = 2, line = 2)
   if (!is.null(fsvsimobj)) {
    if (j <= rtrue) {
     abline(h = fsvsimobj$facload[i,j], col = 3, lty = 2)
    } else {
     abline(h = 0, col = 2, lty = 3)
    }
   }
   if (j %% effrows == 0) {
    axis(1)
    plot.new()
   }
  }
  if (r > effrows) for (k in 1:(effrows - r %% effrows + 1)) if (r %% effrows != 0) plot.new()  # HACK
 }
 par(oldpar)
 invisible(x)
}

#' Density plots of factor loadings draws
#'
#' \code{facloaddensplot} draws kernel smoothed density plots of the marginal
#' factor loadings posterior.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param rows Number of rows per page.
#' @param thesecols Which factor loadings columns should be plotted? Defaults to 1:r.
#' @param xlim Vector of length two containing lower and upper bounds of the
#' horizontal axis. If \code{NULL}, these are automatically determined.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
facloaddensplot <- function(x, fsvsimobj = NULL, rows = 5, thesecols = NULL, xlim = NULL) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (any(dim(x$facload) == 0)) stop("There aren't any factor loadings to plot.")
 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$facload)[2]
 draws <- dim(x$facload)[3]

 if (exists("identifier", x)) {
  ident <- x$identifier[,"identifier"]
 } else {
  ident <- rep(0, r)
 }

 if (r > 1 && all(x$facload[1,2,] == 0)) restrict = "upper" else restrict = "none"

 snames <- colnames(x$y)
 if (is.null(snames)) snames <- 1:m

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
  rtrue <- ncol(fsvsimobj$facload)
 }

 densadjust = 1.7

 if (is.null(thesecols)) thesecols <- 1:r
 if (restrict == "upper") starter <- min(thesecols) else starter <- 1

 oldpar <- par(mfrow = c(rows, length(thesecols)), mgp = c(2, .5, 0), mar = c(1.5, 1.5, 1.5, 0.5))
 for (i in starter:m) {
  for (j in thesecols) {
   if (j > i && restrict == "upper") {
    plot.new()
   } else {
    plot(density(x$facload[i,j,], adjust = densadjust), xlab = "", ylab = "", main = "", xlim = xlim)
    abline(v = 0, lty = 3)
    if (i == ident[j]) cola = 2 else cola = 1
    mtext(paste0("Series ", i, " (", snames[i], ") on Factor ", j), cex = .7, col = cola)
    if (!is.null(fsvsimobj)) {
     if (j <= rtrue) {
      points(fsvsimobj$facload[i,j], 0, col = 3, cex = 1.5, pch = 16)
     } else {
      points(0, 0, col = 2, pch = 18, cex = 1.8)
     }
    }
   }
  }
 }
 par(oldpar)
 invisible(x)
}

#' Several factor SV plots useful for model diagnostics
#'
#' Draws a collection of plots to explore the posterior distribution
#' of a fitted factor SV model.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param ... Other arguments will be passed on to the subfunctions.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
plotalot <- function(x, fsvsimobj = NULL, ...) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 r <- ncol(x$facload)
 if (r > 0) facloadpointplot(x, fsvsimobj = fsvsimobj, ...)
 logvartimeplot(x, fsvsimobj = fsvsimobj, ...)
 if (r > 0) facloaddensplot(x, fsvsimobj = fsvsimobj, ...)
 if (r > 0) facloadtraceplot(x, fsvsimobj = fsvsimobj, ...)
 paratraceplot(x, fsvsimobj = fsvsimobj, ...)
}


#' Default factor SV plot
#'
#' Displays the correlation matrix at the last sampling point in time.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param quantiles Posterior quantiles to be visualized. Must be of length 1 or 3.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an optional object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param col Optional color palette.
#' @param ... Other arguments will be passed on to \link[corrplot]{corrplot}.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
plot.fsvdraws <- function(x, quantiles = c(.05, .5, .95), col = NULL, fsvsimobj = NULL, ...) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
 }
 if (is.null(col)) {
  colpal <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
    "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  col <- rev(colpal(200))
 }
 effdraws <- dim(x$fac)[3]
 n <- nrow(x$y)
 m <- ncol(x$y)

 mycormat <- cormat(x, "last")[,,,1]
 
 # a new version corrplot requires rownames and colnames to be set
 if (is.null(rownames(mycormat))) rownames(mycormat) <- colnames(mycormat) <- seq_len(m)

 if (length(quantiles) == 1) {
   plotCI <- "n"
   mycormatlower <- NULL
   mycormatmed <- apply(mycormat, 1:2, quantile, quantiles)
   mycormatupper <- NULL
 } else if (length(quantiles) == 3) {
   plotCI <- "circle"
   mycormatlower <- apply(mycormat, 1:2, quantile, quantiles[1])
   mycormatmed <- apply(mycormat, 1:2, quantile, quantiles[2])
   mycormatupper <- apply(mycormat, 1:2, quantile, quantiles[3])
 } else {
   stop("Length of argument 'quantiles' must be one or three.")
 }

 corrplot::corrplot(mycormatmed, plotCI = plotCI, lowCI.mat = mycormatlower, uppCI.mat = mycormatupper, col = col)

 if (!is.null(fsvsimobj)) {
   cortrue <- cormat(fsvsimobj, n)[,,1]
   diag(cortrue) <- NA
   oldpar <- par(xpd = TRUE)
   symbols(rep(1:m, each = m), rep(m:1, m),
	   circles = .9*abs(as.numeric(cortrue))^0.5/2,
	   fg = "green", inches = FALSE, add = TRUE)
   par(oldpar)
  }
 invisible(x)
}

#' Plots pairwise correlations over time
#'
#' \code{corplot} gives an overview of (certain) pairwise correlations.
#' Throws a warning if these haven't been stored during sampling.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#' @param fsvsimobj To indicate data generating values in case of simulated
#' data, pass an object of type \code{fsvsim} (usually the result of a
#' call to \code{\link{fsvsim}}).
#' @param these Indicator which correlations should be plotted. Default is all.
#' @param start First point in time to plot.
#' @param end Last point in time to plot.
#' @param maxrows The maximum number of rows per page.
#' @param ... Other arguments will be passed on to \code{\link[stats]{ts.plot}}.
#'
#' @return Returns \code{x} invisibly.
#'
#' @family plotting
#'
#' @export
corplot <- function(x, fsvsimobj = NULL, these = 1:(ncol(x$y)*(ncol(x$y)-1)/2), start = 1,
		    end = nrow(x$y), maxrows = 10, ...) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 if (!exists("runningstore", x) || !exists("cor", x$runningstore)) {
  warning("Correlation hasn't been stored during sampling.")
  return(invisible(x))
 }

 m <- ncol(x$y)
 n <- nrow(x$y)
 r <- dim(x$logvar)[2] - m
 draws <- dim(x$para)[3]
 snames <- colnames(x$y)
 dates <- rownames(x$y)
 if (is.null(dates)) dates <- 1:n

 if (!is.numeric(start) || length(start) != 1 || start > n || start < 1) stop("Illegal argument value 'start'.")
 if (!is.numeric(end) || length(end) != 1 || end > n || end < 1) stop("Illegal argument value 'start'.")
 if (start > end) stop("Argument 'end' must be greater or equal to 'start'.")
 times <- start:end

 if (!is.null(fsvsimobj)) {
  if (!is(fsvsimobj, "fsvsim")) stop("If provided, argument 'fsvsimobj' must be an 'fsvsim' object.")
 }

 if (!is.numeric(these) || min(these) < 1 || max(these) > m*(m-1)/2) stop("Illegal argument value 'these'.")

 oldpar <- par(mfrow = c(min(maxrows, length(these)), 1), mgp = c(2, .5, 0), mar = c(1.5, 1.5, 1.3, 0.5))

 for (i in these) {
  thismean <- x$runningstore$cor[times,i,"mean"]
  thissd <- x$runningstore$cor[times,i,"sd"]
  ts.plot(cbind(thismean - 2*thissd, thismean, thismean + 2*thissd),
	  col = c("gray", 1, "gray"), main = "", xlab = "", ylab = "",
	  gpars = list(xaxt = 'n', ...))
  abline(h = 0, lty = 3)
  ats <- round(seq(1, length(times), len = min(11, length(times))))
  axis(1, labels = dates[times][ats], at = ats)

  whiches <- as.numeric(unlist(strsplit(dimnames(x$runningstore$cor)[[2]][i], "_")))

  title(paste0("Estimated correlation of series ", whiches[2], " (",
	       snames[whiches[2]], ") and series ", whiches[1], " (",
	       snames[whiches[1]], ") (mean +/- 2sd)"))

  if (!is.null(fsvsimobj)) {
   lines(corelement(fsvsimobj, whiches[1], whiches[2], these = times), col = 3)
  }
 }
 par(oldpar)
 invisible(x)
}

#' Plots posterior draws and posterior means of the eigenvalues of crossprod(facload)
#'
#' \code{evdiag} computes, returns, and visualizes the eigenvalues of crossprod(facload).
#' This can be used as a rough guide to choose the numbers of factors in a model.
#'
#' @note Experimental feature. Please be aware that - for the sake of simplicity and
#' interpretability - both the time-varying idiosyncratic as well as the time-varying
#' factor volatilities are simply ignored.
#'
#' @param x Object of class \code{'fsvdraws'}, usually resulting from a call
#' to \code{\link{fsvsample}}.
#'
#' @return Invisibly returns a matrix with posterior samples of the eigenvalues of
#' crossprod(facload)
#'
#' @family plotting
#'
#' @importFrom graphics boxplot
#'
#' @export
evdiag <- function(x) {
 if (!is(x, "fsvdraws")) stop("This function expects an 'fsvdraws' object.")
 l <- x$facload
 d <- matrix(NA_real_, dim(l)[3], dim(l)[2])
 for (i in seq_len(dim(l)[3])) {
  tmp <- l[,,i, drop=FALSE]
  dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2])
  d[i,] <- svd(crossprod(tmp))$d
 }

 xs <- seq_len(dim(l)[2])
 ys <- colMeans(d)
 boxplot(d, border = "grey", main = "Eigenvalues of crossprod(facload)",
	 ylim = c(0, 2*max(ys)))
 lines(xs, ys, lwd = 2, type = "b", pch = 2)
 legend("topright", c("Posterior draws", "Posterior means"), col = c("gray", "black"), lty = 1, pch = c(NA, 2), lwd = c(1, 2))
 invisible(d)
}
