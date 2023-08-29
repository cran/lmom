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
\deqn{F(x)=\exp\lbrace-\exp(-y)\rbrace}{F(x) = exp(-exp(-y))} where
\deqn{y=-k^{-1}\log\lbrace1-k(x-\xi)/\alpha\rbrace,}{y = (-1/k) log(1-k(x-xi)/alpha) ,}
with \eqn{x}{x} bounded by \eqn{\xi+\alpha/k}{xi+alpha/k}
from below if \eqn{k<0}{k<0} and from above if \eqn{k>0}{k>0},
and quantile function
\deqn{x(F)=\xi+{\alpha\over k}\lbrace1-(-\log F)^k\rbrace.}{x(F) = xi + alpha (1 - (-log(F))^k) / k .}

Extreme-value distribution types I, II and III (Gumbel, Frechet, Weibull)
correspond to shape parameter values
\eqn{k=0}{k=0}, \eqn{k<0}{k<0} and \eqn{k>0}{k>0} respectively.
}

\value{\code{cdfgev} gives the distribution function;
       \code{quagev} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.

Two parametrizations of the generalized extreme-value distribution are in common use.
When Jenkinson (1955) introduced the distribution
he wrote the distribution function in the form
\deqn{F(x) = \exp [ - \lbrace 1 - k ( x - \xi ) / \alpha) \rbrace^{1/k}].}{F(x) = exp [ - { 1 - k ( x - xi ) / alpha) } ^ (1/k) ] .}
and that is the form used in \R{} package \pkg{lmom}.  A slight inconvenience with it is that the
skewness of the distribution is a decreasing function of the shape parameter \eqn{k}.
Perhaps for this reason, authors of some other \R{} packages prefer a form in which
the sign of the shape parameter \eqn{k} is changed and the parameters are renamed:
\deqn{F(x) = \exp [ - \lbrace 1 + \xi ( x - \mu ) / \sigma) \rbrace^{-1/\xi}].}{F(x) = exp [ - { 1 + xi ( x - mu ) / sigma) } ^ (-1/xi) ] .}
Users should be able to mix functions from packages that use either form; just be aware that
the sign of the shape parameter will need to be changed when converting from one form to the other
(and that \eqn{\xi}{xi} is a location parameter in one form and a shape parameter in the other).
}

\references{
  Jenkinson, A. F. (1955). The frequency distribution of the annual maximum
  (or minimum) of meteorological elements.
  \emph{Quarterly Journal of the Royal Meteorological Society}, \bold{81}, 158-171.
}



\seealso{
\code{\link{cdfgum}} for the Gumbel (extreme-value type I) distribution.

\code{\link{cdfkap}} for the kappa distribution,
which generalizes the generalized extreme-value distribution.

\code{\link{cdfwei}} for the Weibull distribution,
}

\examples{
# Random sample from the generalized extreme-value distribution
# with parameters xi=0, alpha=1, k=-0.5.
quagev(runif(100), c(0,1,-0.5))
}
\keyword{distribution}
