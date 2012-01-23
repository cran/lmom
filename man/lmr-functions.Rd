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
\tabular{ll}{
\code{lmrexp} \tab exponential \cr
\code{lmrgam} \tab gamma \cr
\code{lmrgev} \tab generalized extreme-value \cr
\code{lmrglo} \tab generalized logistic \cr
\code{lmrgpa} \tab generalized Pareto \cr
\code{lmrgno} \tab generalized normal \cr
\code{lmrgum} \tab Gumbel (extreme-value type I) \cr
\code{lmrkap} \tab kappa \cr
\code{lmrln3} \tab three-parameter lognormal \cr
\code{lmrnor} \tab normal \cr
\code{lmrpe3} \tab Pearson type III \cr
\code{lmrwak} \tab Wakeby \cr
\code{lmrwei} \tab Weibull \cr
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

\author{J. R. M. Hosking \email{hosking@us.ibm.com}}


\seealso{
\code{\link{lmrp}} to compute \eqn{L}-moments of a general distribution
specified by its cumulative distribution function or quantile function.

\code{\link{samlmu}} to compute \eqn{L}-moments of a data sample.

\code{\link{pelexp}}, etc., to compute the parameters
of a distribution given its \eqn{L}-moments.

For individual distributions, see their cumulative distribution functions:
\tabular{ll}{
\code{\link{cdfexp}} \tab exponential \cr
\code{\link{cdfgam}} \tab gamma \cr
\code{\link{cdfgev}} \tab generalized extreme-value \cr
\code{\link{cdfglo}} \tab generalized logistic \cr
\code{\link{cdfgpa}} \tab generalized Pareto \cr
\code{\link{cdfgno}} \tab generalized normal \cr
\code{\link{cdfgum}} \tab Gumbel (extreme-value type I) \cr
\code{\link{cdfkap}} \tab kappa \cr
\code{\link{cdfln3}} \tab three-parameter lognormal \cr
\code{\link{cdfnor}} \tab normal \cr
\code{\link{cdfpe3}} \tab Pearson type III \cr
\code{\link{cdfwak}} \tab Wakeby \cr
\code{\link{cdfwei}} \tab Weibull \cr
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
