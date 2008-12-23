\name{cdfgno}
\alias{cdfgno}
\alias{quagno}
\title{Generalized normal distribution}
\description{Distribution function and quantile function
of the generalized normal distribution.}
\usage{
cdfgno(x, para = c(0, 1, 0))
quagno(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha, k}{xi, alpha, k} (location, scale, shape).}
}

\details{The generalized normal distribution with
location parameter \eqn{\xi}{xi},
scale parameter \eqn{\alpha}{alpha} and
shape parameter \eqn{k}{k} has distribution function
\deqn{F(x)=\Phi(y)}{F(x) = Phi(y)} where
\deqn{y=-k^{-1}\log(1-k(x-\xi)/\alpha)}{y = (-1/k) log(1-k(x-xi)/alpha)}
and \eqn{\Phi(y)}{Phi(y)} is the distribution function of the standard normal distribution,
with \eqn{x}{x} bounded by \eqn{\xi+\alpha/k}{xi+alpha/k}
from below if \eqn{k<0}{k<0} and from above if \eqn{k>0}{k>0}.
%and quantile function
%\deqn{x(F)=\xi-{\alpha\over k}(1-(1-F \over F)^k).}{x(F) = xi - alpha (1 - ((1-F)/F)^k) / k .}

The generalized normal distribution contains as special cases the usual
three-parameter lognormal distribution, corresponding to \eqn{k<0}{k<0},
with a finite lower bound and positive skewness;
the normal distribution, corresponding to \eqn{k=0}{k=0};
and the reverse lognormal distribution, corresponding to \eqn{k>0}{k>0},
with a finite upper bound and negative skewness.
The two-parameter lognormal distribution,
with a lower bound of zero and positive skewness,
is obtained when \eqn{k<0}{k<0} and \eqn{\xi+\alpha/k=0}{xi+alpha/k=0}.
}

\value{\code{cdfgno} gives the distribution function;
       \code{quagno} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfln3}} for the \pkg{lmom} package's version of the three-parameter
lognormal distribution.

\code{\link{cdfnor}} for the \pkg{lmom} package's version of the normal distribution.

\code{\link{pnorm}} for the standard \R version of the normal distribution.

\code{\link{plnorm}} for the standard \R version of the two-parameter lognormal distribution.
}

\examples{
# Random sample from the generalized normal distribution
# with parameters xi=0, alpha=1, k=-0.5.
quagno(runif(100), c(0,1,-0.5))

# The generalized normal distribution with parameters xi=1, alpha=1, k=-1,
# is the standard lognormal distribution.  An illustration:
fval<-seq(0.1,0.9,by=0.1)
cbind(fval, lognormal=qlnorm(fval), g.normal=quagno(fval, c(1,1,-1)))
}
\keyword{distribution}
