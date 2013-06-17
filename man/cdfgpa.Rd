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
distribution functions \code{pnorm}, \code{qnorm}, etc.}


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
