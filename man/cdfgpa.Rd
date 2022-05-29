\name{cdfgpa}
\alias{cdfgpa}
\alias{quagpa}
\title{Generalized Pareto distribution}
\description{Distribution function and quantile function
of the generalized Pareto distribution.}
\usage{
cdfgpa(x, para = c(0, 1, 0))
quagpa(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha, k}{xi, alpha, k} (location, scale, shape).}
}

\details{The generalized Pareto distribution with
location parameter \eqn{\xi}{xi},
scale parameter \eqn{\alpha}{alpha} and
shape parameter \eqn{k}{k} has distribution function
\deqn{F(x)=1-\exp(-y)}{F(x) = 1 - exp(-y)} where
\deqn{y=-k^{-1}\log\lbrace1-k(x-\xi)/\alpha\rbrace,}{y = (-1/k) log(1-k(x-xi)/alpha) ,}
with \eqn{x}{x} bounded by \eqn{\xi+\alpha/k}{xi+alpha/k}
from below if \eqn{k<0}{k<0} and from above if \eqn{k>0}{k>0},
and quantile function
\deqn{x(F)=\xi+{\alpha\over k}\lbrace 1-(1-F)^k\rbrace.}{x(F) = xi + alpha (1 - (1-F)^k) / k .}

The exponential distribution is the special case \eqn{k=0}{k=0}.
The uniform distribution is the special case \eqn{k=1}{k=1}.
}

\value{\code{cdfgpa} gives the distribution function;
       \code{quagpa} gives the quantile function.}

\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.

Two parametrizations of the generalized Pareto distribution are in common use.
When Jenkinson (1955) introduced the generalized extreme-value distribution
he wrote the distribution function in the form
\deqn{F(x) = \exp [ - \lbrace 1 - k ( x - \xi ) / \alpha) \rbrace^{1/k}].}{F(x) = exp [ - { 1 - k ( x - xi ) / alpha) } ^ (1/k) ] .}
Hosking and Wallis (1987) wrote the distribution function of the
generalized Pareto distribution analogously as
\deqn{F(x) = 1 - \lbrace 1 - k ( x - \xi ) / \alpha) \rbrace^{1/k}}{F(x) = 1 - { 1 - k ( x - xi ) / alpha) } ^ (1/k)}
and that is the form used in \R{} package \pkg{lmom}.  A slight inconvenience with it is that the
skewness of the distribution is a decreasing function of the shape parameter \eqn{k}.
Perhaps for this reason, authors of some other \R{} packages prefer a form in which
the sign of the shape parameter \eqn{k} is changed and the parameters are renamed:
\deqn{F(x) = 1 - \lbrace 1 + \xi ( x - \mu ) / \sigma) \rbrace ^{-1/\xi}.}{F(x) = 1 - { 1 - xi ( x - mu ) / sigma) } ^ (-1/xi) .}
Users should be able to mix functions from packages that use either form; just be aware that
the sign of the shape parameter will need to be changed when converting from one form to the other
(and that \eqn{\xi}{xi} is a location parameter in one form and a shape parameter in the other).
}

\references{
  Hosking, J. R. M., and Wallis, J. R. (1987). Parameter and quantile estimation
  for the generalized Pareto distribution.
  \emph{Technometrics}, \bold{29}, 339-349.

  Jenkinson, A. F. (1955). The frequency distribution of the annual maximum
  (or minimum) of meteorological elements.
  \emph{Quarterly Journal of the Royal Meteorological Society}, \bold{81}, 158-171.
}

\seealso{
\code{\link{cdfexp}} for the exponential distribution.

\code{\link{cdfkap}} for the kappa distribution and
\code{\link{cdfwak}} for the Wakeby distribution,
which generalize the generalized Pareto distribution.
}

\examples{
# Random sample from the generalized Pareto distribution
# with parameters xi=0, alpha=1, k=-0.5.
quagpa(runif(100), c(0,1,-0.5))
}
\keyword{distribution}
