\name{cdfglo}
\alias{cdfglo}
\alias{quaglo}
\title{Generalized logistic distribution}
\description{Distribution function and quantile function
of the generalized logistic distribution.}
\usage{
cdfglo(x, para = c(0, 1, 0))
quaglo(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha, k}{xi, alpha, k} (location, scale, shape).}
}

\details{The generalized logistic distribution with
location parameter \eqn{\xi}{xi},
scale parameter \eqn{\alpha}{alpha} and
shape parameter \eqn{k}{k} has distribution function
\deqn{F(x)=1/\lbrace 1+\exp(-y)\rbrace}{F(x)= 1 / (1+exp(-y))} where
\deqn{y=-k^{-1}\log\lbrace1-k(x-\xi)/\alpha\rbrace,}{y = (-1/k) log(1-k(x-xi)/alpha) ,}
with \eqn{x}{x} bounded by \eqn{\xi+\alpha/k}{xi+alpha/k}
from below if \eqn{k<0}{k<0} and from above if \eqn{k>0}{k>0},
and quantile function
%\deqn{x(F)=\xi+{\alpha\over k}(1-({1-F \over F})^k).}{x(F) = xi + alpha (1 - ((1-F)/F)^k) / k .}
\deqn{x(F)=\xi+{\alpha\over k}\biggl\lbrace1-\biggl({1-F \over F}\biggr)^k\biggr\rbrace.}{x(F) = xi + alpha (1 - ((1-F)/F)^k) / k .}

The logistic distribution is the special case \eqn{k=0}{k=0}.
}

\value{\code{cdfglo} gives the distribution function;
       \code{quaglo} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfkap}} for the kappa distribution,
which generalizes the generalized logistic distribution.
}

\examples{
# Random sample from the generalized logistic distribution
# with parameters xi=0, alpha=1, k=-0.5.
quaglo(runif(100), c(0,1,-0.5))
}
\keyword{distribution}
