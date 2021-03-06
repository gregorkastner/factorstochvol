% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{evdiag}
\alias{evdiag}
\title{Plots posterior draws and posterior means of the eigenvalues of crossprod(facload)}
\usage{
evdiag(x)
}
\arguments{
\item{x}{Object of class \code{'fsvdraws'}, usually resulting from a call
to \code{\link{fsvsample}}.}
}
\value{
Invisibly returns a matrix with posterior samples of the eigenvalues of
crossprod(facload)
}
\description{
\code{evdiag} computes, returns, and visualizes the eigenvalues of crossprod(facload).
This can be used as a rough guide to choose the numbers of factors in a model.
}
\note{
Experimental feature. Please be aware that - for the sake of simplicity and
interpretability - both the time-varying idiosyncratic as well as the time-varying
factor volatilities are simply ignored.
}
\seealso{
Other plotting: 
\code{\link{comtimeplot}()},
\code{\link{corimageplot}()},
\code{\link{corplot}()},
\code{\link{cortimeplot}()},
\code{\link{facloadcredplot}()},
\code{\link{facloaddensplot}()},
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
