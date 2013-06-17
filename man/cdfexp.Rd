\name{cdfexp}
\alias{cdfexp}
\alias{quaexp}
\title{Exponential distribution}
\description{Distribution function and quantile function
of the exponential distribution.}
\usage{
cdfexp(x, para = c(0, 1))
quaexp(f, para = c(0, 1))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha}{xi, alpha} (location, scale).}
}

\details{The exponential distribution with parameters
\eqn{\xi}{xi} (lower bound) and \eqn{\alpha}{alpha} (mean)
has distribution function
\deqn{F(x)=1-\exp\lbrace-(x-\xi)/\alpha\rbrace}{F(x) = 1 - exp(-(x-xi)/alpha)}
for \eqn{x\ge0}{x >= 0}, and quantile function
\deqn{x(F)=\xi-\alpha\log(1-F).}{x(F) = xi - alpha log(1-F) .}
}

\value{\code{cdfexp} gives the distribution function;
       \code{quaexp} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R functions
\code{pexp} and \code{qexp}.
}


\seealso{
\code{\link{pexp}} for the standard \R version of the exponential distribution.
\code{\link{cdfgam}} for the gamma distribution,
\code{\link{cdfgpa}} for the generalized Pareto distribution,
\code{\link{cdfkap}} for the kappa distribution,
\code{\link{cdfpe3}} for the Pearson type III distribution,
and \code{\link{cdfwak}} for the Wakeby distribution,
all of which generalize the exponential distribution.
}

\examples{
# Random sample from the exponential distribution
# with lower bound 0 and mean 3.
quaexp(runif(100), c(0,3))
}
\keyword{distribution}
