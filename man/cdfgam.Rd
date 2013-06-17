\name{cdfgam}
\alias{cdfgam}
\alias{quagam}
\title{Gamma distribution}
\description{Distribution function and quantile function
of the gamma distribution.}
\usage{
cdfgam(x, para = c(1, 1))
quagam(f, para = c(1, 1))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\alpha, \beta}{alpha, beta} (shape, scale).}
}

\details{The gamma distribution with
shape parameter \eqn{\alpha}{alpha} and
scale parameter \eqn{\beta}{beta} has probability density function
\deqn{f(x)={x^{\alpha-1} \exp(-x/\beta) \over \beta^\alpha \Gamma(\alpha)}}{f(x) = x^alpha exp(-x/beta) / ( beta^alpha Gamma(alpha) )}
for \eqn{x\ge0}{x >= 0}, where \eqn{\Gamma(.)}{Gamma(.)} is the gamma function.
}

\value{\code{cdfgam} gives the distribution function;
       \code{quagam} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R functions
\code{pgamma} and \code{qgamma}.
}


\seealso{
\code{\link{gamma}} for the gamma function.

\code{\link{pgamma}} for the standard \R version of the gamma distribution.

\code{\link{cdfpe3}} for the Pearson type III distribution,
which generalizes the gamma distribution.
}

\examples{
# Random sample from the gamma distribution
# with shape parameter 4 and mean 1.
quagam(runif(100), c(4,1/4))
}
\keyword{distribution}
