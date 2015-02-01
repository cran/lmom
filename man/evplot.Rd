\name{evplot}
\alias{evplot}
\alias{evplot.default}
\alias{evdistp}
\alias{evdistq}
\alias{evpoints}
\title{Extreme-value plot}
\description{
\code{evplot} draws an \dQuote{extreme-value plot}, i.e. a quantile-quantile plot
in which the horizontal axis is the quantile of an
extreme-value type I (Gumbel) distribution.

\code{evdistp} adds the cumulative distribution function of a distribution
to an extreme-value plot.

\code{evdistq} adds the quantile function of a distribution
to an extreme-value plot.

\code{evpoints} adds a set of data points to an extreme-value plot.
}
\usage{
evplot(y, ...)

\method{evplot}{default}(y, qfunc, para, npoints = 101, plim, xlim = c(-2, 5),
       ylim, type,
       xlab = expression("Reduced variate,  " * -log(-log(italic(F)))),
       ylab = "Quantile", rp.axis = TRUE, ...)

evdistp(pfunc, para, npoints = 101, ...)
evdistq(qfunc, para, npoints = 101, ...)

evpoints(y, ...)
}
\arguments{
\item{y}{Numeric vector.  The data values in the vector are plotted
  on the extreme-value plot.}
\item{qfunc}{A quantile function.  The function is drawn as a curve
  on the extreme-value plot.}
\item{pfunc}{A cumulative distribution function.  The function is drawn as a curve
  on the extreme-value plot.}
\item{para}{Distribution parameters for the quantile function \code{qfunc}
  or cumulative distribution function \code{pfunc}.

  If \code{pfunc} or \code{qfunc} is the standard \R form of quantile function,
  \code{para} should be a list.

  If \code{pfunc} or \code{qfunc} is the \code{qua...} form
  of quantile function used throughout the \pkg{lmom} package,
  \code{para} should be a numeric vector.

  In \code{evplot}, \code{para} is not used if \code{qfunc} is omitted.}
\item{npoints}{Number of points to use in drawing the quantile function.
  The points are equally spaced along the x axis.
  Not used if \code{qfunc} is omitted.}
\item{plim}{X axis limits, specified as probabilities.}
\item{xlim}{X axis limits, specified as values of the Gumbel reduced variate
  \eqn{-\log(-\log F)}{-log(-log(F))}, where \eqn{F}{F} is the
  nonexceedance probability.
  Not used if \code{plim} is specified.}
\item{ylim}{Y axis limits.}
\item{type}{Plot type.  Determines how the data values in \code{y} are plotted.
  Interpreted in the same way as argument \code{type} of function \code{plot},
  i.e. \code{"p"} for points, \code{"b"} for points connected by lines, etc.}
\item{xlab}{X axis label.}
\item{ylab}{Y axis label.}
\item{rp.axis}{Logical. Whether to draw the \dQuote{Return period} axis,
  a secondary horizontal axis.}
\item{\dots}{Additional arguments are passed to the plotting routine.}
}
%\details{ }
%\value{ }

\note{Data points are plotted at the Gringorten plotting position,
i.e. the \eqn{i}{i'}th smallest of \eqn{n}{n} data points is plotted
at the horizontal position corresponding to nonexceedance probability
\eqn{(i-0.44)/(n+0.12)}.
}

\section{Arguments of cumulative distribution functions and quantile functions}{
\code{pfunc} and \code{qfunc} can be either the standard \R form of
cumulative distribution function or quantile function
(i.e. for a distribution with \eqn{r} parameters, the first argument is the
variate \eqn{x} or the probability \eqn{p} and the next \eqn{r} arguments
are the parameters of the distribution) or the \code{cdf\dots} or
\code{qua\dots} forms used throughout the \pkg{lmom} package
(i.e. the first argument is the variate \eqn{x} or probability \eqn{p}
and the second argument is a vector containing the parameter values).
}

\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}

\examples{
# Extreme-value plot of Ozone from the airquality data
data(airquality)
evplot(airquality$Ozone)

# Fit a GEV distribution and add it to the plot
evdistq(quagev, pelgev(samlmu(airquality$Ozone)))

# Not too good -- try a kappa distribution instead
evdistq(quakap, pelkap(samlmu(airquality$Ozone)), col="red")
}
\keyword{hplot}
