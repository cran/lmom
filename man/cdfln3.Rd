\name{cdfln3}
\alias{cdfln3}
\alias{qualn3}
\title{Three-parameter lognormal distribution}
\description{Distribution function and quantile function
of the three-parameter lognormal distribution.}
\usage{
cdfln3(x, para = c(0, 0, 1))
qualn3(f, para = c(0, 0, 1))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\zeta, \mu, \sigma}{zeta, mu, sigma}
(lower bound, mean on log scale, standard deviation on log scale).}
}

\details{The three-parameter lognormal distribution with
lower bound \eqn{\zeta}{zeta},
mean on log scale \eqn{\mu}{mu}, and
standard deviation on log scale \eqn{\sigma}{sigma} has distribution function
\deqn{F(x)=\Phi(y),}{F(x) = Phi(y),} \eqn{x>0}, where
\deqn{y=\lbrace\log(x - \zeta)-\mu\rbrace/\sigma}{y = (log(x-zeta) - mu) / sigma}
and \eqn{\Phi(y)}{Phi(y)} is the distribution function of the standard
normal distribution.
}

\value{\code{cdfln3} gives the distribution function;
       \code{qualn3} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgno}} for the generalized normal distribution,
a more general form of the three-parameter lognormal distribution.

\code{\link{qlnorm}} for the standard \R version of the
two-parameter lognormal distribution.
}

\examples{
# Random sample from three-parameter lognormal distribution
# with parameters zeta=0, mu=1, sigma=0.5.
qualn3(runif(100), c(0,1,0.5))

## Functions for the three-parameter lognormal distribution can
## also be used with the two-parameter lognormal distribution
# Generate a random sample from a standard lognormal distribution
xx <- qualn3(runif(50))
# Fit 2-parameter LN distribution
pelln3(samlmu(xx), bound=0)
# Fit 2-parameter LN distribution "in log space",
# i.e. fit normal distribution to log-transformed data
pelnor(samlmu(log(xx)))
}
\keyword{distribution}
