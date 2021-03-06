% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{facloaddensplot}
\alias{facloaddensplot}
\title{Density plots of factor loadings draws}
\usage{
facloaddensplot(x, fsvsimobj = NULL, rows = 5, thesecols = NULL, xlim = NULL)
}
\arguments{
\item{x}{Object of class \code{'fsvdraws'}, usually resulting from a call
to \code{\link{fsvsample}}.}

\item{fsvsimobj}{To indicate data generating values in case of simulated
data, pass an object of type \code{fsvsim} (usually the result of a
call to \code{\link{fsvsim}}).}

\item{rows}{Number of rows per page.}

\item{thesecols}{Which factor loadings columns should be plotted? Defaults to 1:r.}

\item{xlim}{Vector of length two containing lower and upper bounds of the
horizontal axis. If \code{NULL}, these are automatically determined.}
}
\value{
Returns \code{x} invisibly.
}
\description{
\code{facloaddensplot} draws kernel smoothed density plots of the marginal
factor loadings posterior.
}
\seealso{
Other plotting: 
\code{\link{comtimeplot}()},
\code{\link{corimageplot}()},
\code{\link{corplot}()},
\code{\link{cortimeplot}()},
\code{\link{evdiag}()},
\code{\link{facloadcredplot}()},
\code{\link{facloadpairplot}()},
\code{\link{facloadpointplot}()},
\code{\link{facloadtraceplot}()},
\code{\link{logvartimeplot}()},
\code{\link{paratraceplot}()},
\code{\link{plot.fsvdraws}()},
\code{\link{plotalot}()},
\code{\link{voltimeplot}()}
}
\concept{plotting}
