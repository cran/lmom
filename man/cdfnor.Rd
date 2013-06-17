\name{cdfnor}
\alias{cdfnor}
\alias{quanor}
\title{Normal distribution}
\description{Distribution function and quantile function
of the normal distribution.}
\usage{
cdfnor(x, para = c(0, 1))
quanor(f, para = c(0, 1))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\mu, \sigma}{mu, sigma} (location, scale).}
}

\details{The normal distribution with
location parameter \eqn{\mu}{mu} and
scale parameter \eqn{\sigma}{sigma}
has probability density function
\deqn{f(x)={1\over \sigma\sqrt{2\pi}} \exp\lbrace-(x-\mu)^2/(2 \sigma^2)\rbrace.}{f(x) = (1/ sigma sqrt(2 pi)) exp( -(x-mu)^2 / (2 sigma^2) ) .}
}

\value{\code{cdfnor} gives the distribution function;
       \code{quanor} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm} and \code{qnorm}.}


\seealso{
\code{\link{pnorm}} for the standard \R version of the normal distribution.

\code{\link{cdfgno}} for the generalized normal distribution,
which generalizes the normal distribution.
}

\examples{
# Random sample from the normal distribution
# with mean 0 and standard deviation 3.
quanor(runif(100), c(0,3))
}
\keyword{distribution}
