\name{samlmu}
\alias{samlmu}
\alias{samlmu.s}
\title{Sample L-moments}
\description{
Computes the \dQuote{unbiased} sample \eqn{L}-moments and \eqn{L}-moment ratios of a data vector.
}
\usage{
samlmu(x, nmom = 4, sort.data = TRUE)
samlmu.s(x, nmom = 4, sort.data = TRUE)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{nmom}{Number of \eqn{L}-moments to be found.}
  \item{sort.data}{Whether the \code{x} vector needs to be sorted.}
}
\details{\code{samlmu} and \code{samlmu.s} are functionally identical.
\code{samlmu} calls a Fortran routine internally, and is usually faster.
\code{samlmu.s} is written entirely in the S language; it is provided
so that users can conveniently see how the calculations are done.}
\value{
  A numeric vector containing the \eqn{L}-moments and \eqn{L}-moment ratios,
  in the order \eqn{\ell_1}{l_1}, \eqn{\ell_2}{l_2},
  \eqn{t_3}{t_3}, \eqn{t_4}{t_4}, etc.
}

\note{Sample \eqn{L}-moments are defined in Hosking (1990).
Calculations use the algorithm given in Hosking (1996, p.14).}

\references{
  Hosking, J. R. M. (1990).
  \eqn{L}-moments: analysis and estimation of distributions
  using linear combinations of order statistics.
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{52}, 105-124.

  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}

\author{J. R. M. Hosking \email{hosking@us.ibm.com}}

\examples{
    data(airquality)
    samlmu(airquality$Ozone,6)
}
\keyword{univar}
