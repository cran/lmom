\name{cdfgev}
\alias{cdfgev}
\alias{quagev}
\title{Generalized extreme-value distribution}
\description{Distribution function and quantile function
of the generalized extreme-value distribution.}
\usage{
cdfgev(x, para = c(0, 1, 0))
quagev(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha, k}{xi, alpha, k} (location, scale, shape).}
}

\details{The generalized extreme-value distribution with
location parameter \eqn{\xi}{xi},
scale parameter \eqn{\alpha}{alpha} and
shape parameter \eqn{k}{k} has distribution function
\deqn{F(x)=\exp(-\exp(-y))}{F(x) = exp(-exp(-y))} where
\deqn{y=-k^{-1}\log(1-k(x-\xi)/\alpha),}{y = (-1/k) log(1-k(x-xi)/alpha) ,}
with \eqn{x}{x} bounded by \eqn{\xi+\alpha/k}{xi+alpha/k}
from below if \eqn{k<0}{k<0} and from above if \eqn{k>0}{k>0},
and quantile function
\deqn{x(F)=\xi-{\alpha\over k}(1-(-\log F)^k).}{x(F) = xi - alpha (1 - (-log(F))^k) / k .}

Extreme-value distribution types I, II and III (Gumbel, Frechet, Weibull)
correspond to shape parameter values
\eqn{k=0}{k=0}, \eqn{k<0}{k<0} and \eqn{k>0}{k>0} respectively.
}

\value{\code{cdfgev} gives the distribution function;
       \code{quagev} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgum}} for the Gumbel (extreme-value type I) distribution.

\code{\link{cdfkap}} for the kappa distribution,
which generalizes the generalized extreme-value distribution.
}

\examples{
# Random sample from the generalized extreme-value distribution
# with parameters xi=0, alpha=1, k=-0.5.
quagev(runif(100), c(0,1,-0.5))
}
\keyword{distribution}
