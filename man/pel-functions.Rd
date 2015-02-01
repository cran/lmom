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
\tabular{rll}{
\tab \code{pelexp} \tab exponential \cr
\tab \code{pelgam} \tab gamma \cr
\tab \code{pelgev} \tab generalized extreme-value \cr
\tab \code{pelglo} \tab generalized logistic \cr
\tab \code{pelgpa} \tab generalized Pareto \cr
\tab \code{pelgno} \tab generalized normal \cr
\tab \code{pelgum} \tab Gumbel (extreme-value type I) \cr
\tab \code{pelkap} \tab kappa \cr
\tab \code{pelln3} \tab three-parameter lognormal \cr
\tab \code{pelnor} \tab normal \cr
\tab \code{pelpe3} \tab Pearson type III \cr
\tab \code{pelwak} \tab Wakeby \cr
\tab \code{pelwei} \tab Weibull \cr
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
Hosking (1996, pp. 10--11).
Exception:
if \code{pelwak} is unable to fit a Wakeby distribution using all 5 \eqn{L}-moments,
it instead fits a generalized Pareto distribution to the first 3 \eqn{L}-moments.
(The corresponding routine in the LMOMENTS Fortran package
would attempt to fit a Wakeby distribution with lower bound zero.)
}

\value{A numeric vector containing the parameters of the distribution.}

\references{
  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}

\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}

\seealso{
\code{\link{pelp}} for parameter estimation of a general distribution
specified by its cumulative distribution function or quantile function.

\code{\link{lmrexp}}, etc., to compute the \eqn{L}-moments
of a distribution given its parameters.

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
# Sample L-moments of Ozone from the airquality data
data(airquality)
lmom <- samlmu(airquality$Ozone)

# Fit a GEV distribution
pelgev(lmom)
}
\keyword{models}
