\name{pel-functions}
\alias{pel...}
\alias{pelexp}
\alias{pelgam}
\alias{pelgev}
\alias{pelglo}
\alias{pelgno}
\alias{pelgpa}
\alias{pelgum}
\alias{pelkap}
\alias{pelln3}
\alias{pelnor}
\alias{pelpe3}
\alias{pelwak}
\alias{pelwei}
\title{Parameter estimation for specific distributions by the method of L-moments}
\description{Computes the parameters of a probability distribution
as a function of the \eqn{L}-moments.
The following distributions are recognized:
\tabular{ll}{
\code{pelexp} \tab exponential \cr
\code{pelgam} \tab gamma \cr
\code{pelgev} \tab generalized extreme-value \cr
\code{pelglo} \tab generalized logistic \cr
\code{pelgpa} \tab generalized Pareto \cr
\code{pelgno} \tab generalized normal \cr
\code{pelgum} \tab Gumbel (extreme-value type I) \cr
\code{pelkap} \tab kappa \cr
\code{pelln3} \tab three-parameter lognormal \cr
\code{pelnor} \tab normal \cr
\code{pelpe3} \tab Pearson type III \cr
\code{pelwak} \tab Wakeby \cr
\code{pelwei} \tab Weibull \cr
}
}
\usage{
pelexp(lmom)
pelgam(lmom)
pelgev(lmom)
pelglo(lmom)
pelgno(lmom)
pelgpa(lmom, bound = NULL)
pelgum(lmom)
pelkap(lmom)
pelln3(lmom, bound = NULL)
pelnor(lmom)
pelpe3(lmom)
pelwak(lmom, bound = NULL, verbose = FALSE)
pelwei(lmom, bound = NULL)
}

\arguments{
\item{lmom}{Numeric vector containing the \eqn{L}-moments of the distribution
or of a data sample.}
\item{bound}{Lower bound of the distribution.  If \code{NULL} (the default),
             the lower bound will be estimated along with the other parameters.}
\item{verbose}{Logical: whether to print a message when not all
               parameters of the distribution can be computed.}
}

\details{Numerical methods and accuracy are as described in
Hosking (1996, pp. 8--9).}

\value{A numeric vector containing the parameters of the distribution.}

\references{
  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}

\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}

\seealso{
\code{\link{pelp}} for parameter estimation of a general distribution
specified by its cumulative distribution function or quantile function.

\code{\link{lmrexp}}, etc., to compute the \eqn{L}-moments
of a distribution given its parameters.

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
# Sample L-moments of Ozone from the airquality data
data(airquality)
lmom <- samlmu(airquality$Ozone)

# Fit a GEV distribution
pelgev(lmom)
}
\keyword{models}
