\name{lmr-functions}
\alias{lmr...}
\alias{lmrexp}
\alias{lmrgam}
\alias{lmrgev}
\alias{lmrglo}
\alias{lmrgno}
\alias{lmrgpa}
\alias{lmrgum}
\alias{lmrkap}
\alias{lmrln3}
\alias{lmrnor}
\alias{lmrpe3}
\alias{lmrwak}
\alias{lmrwei}
\title{L-moments of specific probability distributions}
\description{Computes the \eqn{L}-moments of a probability distribution
given its parameters.
The following distributions are recognized:
\tabular{rll}{
\tab \code{lmrexp} \tab exponential \cr
\tab \code{lmrgam} \tab gamma \cr
\tab \code{lmrgev} \tab generalized extreme-value \cr
\tab \code{lmrglo} \tab generalized logistic \cr
\tab \code{lmrgpa} \tab generalized Pareto \cr
\tab \code{lmrgno} \tab generalized normal \cr
\tab \code{lmrgum} \tab Gumbel (extreme-value type I) \cr
\tab \code{lmrkap} \tab kappa \cr
\tab \code{lmrln3} \tab three-parameter lognormal \cr
\tab \code{lmrnor} \tab normal \cr
\tab \code{lmrpe3} \tab Pearson type III \cr
\tab \code{lmrwak} \tab Wakeby \cr
\tab \code{lmrwei} \tab Weibull \cr
}
}
\usage{
lmrexp(para = c(0, 1), nmom = 2)
lmrgam(para = c(1, 1), nmom = 2)
lmrgev(para = c(0, 1, 0), nmom = 3)
lmrglo(para = c(0, 1, 0), nmom = 3)
lmrgno(para = c(0, 1, 0), nmom = 3)
lmrgpa(para = c(0, 1, 0), nmom = 3)
lmrgum(para = c(0, 1), nmom = 2)
lmrkap(para = c(0, 1, 0, 0), nmom = 4)
lmrln3(para = c(0, 0, 1), nmom = 3)
lmrnor(para = c(0, 1), nmom = 2)
lmrpe3(para = c(0, 1, 0), nmom = 3)
lmrwak(para = c(0, 1, 0, 0, 0), nmom = 5)
lmrwei(para = c(0, 1, 1), nmom = 3)
}

\arguments{
\item{para}{Numeric vector containing the parameters of the distribution.}
\item{nmom}{The number of \eqn{L}-moments to be calculated.}
}

\details{Numerical methods and accuracy are as described in
Hosking (1996, pp. 8--9).}

\value{Numeric vector containing the \eqn{L}-moments.}

\references{
  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}

\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}

\seealso{
\code{\link{lmrp}} to compute \eqn{L}-moments of a general distribution
specified by its cumulative distribution function or quantile function.

\code{\link{samlmu}} to compute \eqn{L}-moments of a data sample.

\code{\link{pelexp}}, etc., to compute the parameters
of a distribution given its \eqn{L}-moments.

For individual distributions, see their cumulative distribution functions:
\tabular{rll}{
\tab \code{\link{cdfexp}} \tab exponential \cr
\tab \code{\link{cdfgam}} \tab gamma \cr
\tab \code{\link{cdfgev}} \tab generalized extreme-value \cr
\tab \code{\link{cdfglo}} \tab generalized logistic \cr
\tab \code{\link{cdfgpa}} \tab generalized Pareto \cr
\tab \code{\link{cdfgno}} \tab generalized normal \cr
\tab \code{\link{cdfgum}} \tab Gumbel (extreme-value type I) \cr
\tab \code{\link{cdfkap}} \tab kappa \cr
\tab \code{\link{cdfln3}} \tab three-parameter lognormal \cr
\tab \code{\link{cdfnor}} \tab normal \cr
\tab \code{\link{cdfpe3}} \tab Pearson type III \cr
\tab \code{\link{cdfwak}} \tab Wakeby \cr
\tab \code{\link{cdfwei}} \tab Weibull \cr
}
}

\examples{
# Compare sample L-moments of Ozone from the airquality data
# with the L-moments of a GEV distribution fitted to the data
data(airquality)
smom <- samlmu(airquality$Ozone, nmom=6)
gevpar <- pelgev(smom)
pmom <- lmrgev(gevpar, nmom=6)
print(smom)
print(pmom)
}
\keyword{distribution}
