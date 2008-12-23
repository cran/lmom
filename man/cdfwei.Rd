\name{cdfwei}
\alias{cdfwei}
\alias{quawei}
\title{Weibull distribution}
\description{Distribution function and quantile function
of the Weibull distribution.}
\usage{
cdfwei(x, para = c(0, 1, 0))
quawei(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\zeta, \beta, \delta}{zeta, beta, delta}
(location, scale, shape).}
}

\details{The Weibull distribution with
location parameter \eqn{\zeta}{zeta},
scale parameter \eqn{\beta}{beta} and
shape parameter \eqn{\delta}{\delta} has distribution function
\deqn{F(x)=1-\exp[-\lbrace(x-\zeta)/\beta\rbrace^\delta]}{F(x) = 1 - exp[ - \{ (x - zeta) /beta \}^delta ]}
for \eqn{x>\zeta}{x>zeta}.
}

\value{\code{cdfwei} gives the distribution function;
       \code{quawei} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgev}} for the generalized extreme-value distribution,
of which the Weibull (reflected through the origin) is a special case.
}

\examples{
# Random sample from a 2-parameter Weibull distribution
# with scale parameter 2 and shape parameter 1.5.
quagev(runif(100), c(0,2,-0.5))

# Illustrate the relation between Weibull and GEV distributions.
# weifit() fits a Weibull distribution to data and returns
#   quantiles of the fitted distribution
# gevfit() fits a Weibull distribution as a "reverse GEV",
#   i.e. fits a GEV distribution to the negated data,
#   then computes negated quantiles
weifit <- function(qval, x) quawei(qval, pelwei(samlmu(x)))
gevfit <- function(qval, x) -quagev(1-qval, pelgev(samlmu(-x)))
# Compare on Ozone data
data(airquality)
weifit(c(0.2,0.5,0.8), airquality$Ozone)
gevfit(c(0.2,0.5,0.8), airquality$Ozone)
}
\keyword{distribution}
