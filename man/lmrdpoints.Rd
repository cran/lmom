\name{lmrdpoints}
\alias{lmrdpoints}
\alias{lmrdlines}
\title{Add points or lines to an L-moment ratio diagram}
\description{
\code{lmrdpoints} adds points,
and \code{lmrdlines} adds connected line segments,
to an \eqn{L}-moment ratio diagram.
}
\usage{
lmrdpoints(x, y=NULL, type="p", ...)
lmrdlines(x, y=NULL, type="l", ...)
}
\arguments{
  \item{x}{Numeric vector of \eqn{L}-skewness values.}
  \item{y}{Numeric vector of \eqn{L}-kurtosis values.
           May be omitted: see \dQuote{Details} below.}
  \item{type}{Character indicating the type of plotting.
              Can be any valid value for the \code{type} argument
              of \code{plot.default}.}
  \item{...}{Further arguments (graphics parameters),
             passed to \code{points} or \code{lines}.}
}
\details{
The functions \code{lmrdpoints} and \code{lmrdlines} are equivalent to
\code{points} and \code{lines} respectively,
except that if argument \code{y} is omitted, \code{x} is assumed to be
an object that contains both \eqn{L}-skewness and \eqn{L}-kurtosis values.
As in \code{lmrd}, it can be a vector with elements named
\code{"t_3"} and \code{"t_4"} (or \code{"tau_3"} and \code{"tau_4"}),
a matrix or data frame with columns named
\code{"t_3"} and \code{"t_4"} (or \code{"tau_3"} and \code{"tau_4"}),
or an object of class \code{"regdata"} (as defined in package \pkg{lmomRFA}).
}
%\value{}
%\references{}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
%\note{}
\seealso{\code{\link{lmrd}}, \code{\link{points}}, \code{\link{lines}}.}
\examples{
# Plot L-moment ratio diagram of Wind from the airquality data set
data(airquality)
lmrd(samlmu(airquality$Wind), xlim=c(-0.2, 0.2))
# Sample L-moments of each month's data
( lmom.monthly <- with(airquality,
  t(sapply(5:9, function(mo) samlmu(Wind[Month==mo])))) )
# Add the monthly values to the plot
lmrdpoints(lmom.monthly, pch=19, col="blue")


# Draw an L-moment ratio diagram and add a line for the
# Weibull distribution
lmrd(xaxs="i", yaxs="i", las=1)
weimom <- sapply( seq(0, 0.9, by=0.01),
  function(tau3) lmrwei(pelwei(c(0,1,tau3)), nmom=4) )
lmrdlines(t(weimom), col='darkgreen', lwd=2)

}
\keyword{hplot}
