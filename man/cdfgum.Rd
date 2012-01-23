\name{cdfgum}
\alias{cdfgum}
\alias{quagum}
\title{Gumbel (extreme-value type I) distribution}
\description{Distribution function and quantile function
of the Gumbel distribution.}
\usage{
cdfgum(x, para = c(0, 1))
quagum(f, para = c(0, 1))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha}{xi, alpha} (location, scale).}
}

\details{The Gumbel distribution with
location parameter \eqn{\xi}{xi} and
scale parameter \eqn{\alpha}{alpha}
has distribution function
\deqn{F(x)=\exp[-\exp\lbrace-(x-\xi)/\alpha\rbrace]}{F(x) = exp(-exp(-(x-xi)/alpha))}
and quantile function
\deqn{x(F)=\xi-\alpha\log(-\log F).}{x(F) = xi - alpha log(-log F) .}
}

\value{\code{cdfgum} gives the distribution function;
       \code{quagum} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgev}} for the generalized extreme-value distribution,
which generalizes the Gumbel distribution.
}

\examples{
# Random sample from the Gumbel distribution with parameters xi=0, alpha=3.
quagum(runif(100), c(0,3))
}
\keyword{distribution}
