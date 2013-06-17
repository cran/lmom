\name{samlmu}
\alias{samlmu}
\alias{samlmu.s}
\alias{.samlmu}
\title{Sample L-moments}
\description{
Computes the \dQuote{unbiased} sample (trimmed) \eqn{L}-moments
and \eqn{L}-moment ratios of a data vector.
}
\usage{
samlmu(x, nmom=4, sort.data=TRUE, ratios=sort.data, trim=0)
samlmu.s(x, nmom=4, sort.data=TRUE, ratios=sort.data, trim=0)
.samlmu(x, nmom=4)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{nmom}{Number of \eqn{L}-moments to be found.}
  \item{sort.data}{Logical: whether the \code{x} vector needs to be sorted.}
  \item{ratios}{Logical. If \code{FALSE}, \eqn{L}-moments are computed;
    if \code{TRUE} (the default), \eqn{L}-moment ratios are computed.}
  \item{trim}{Degree of trimming.
    If a single value, symmetric trimming of the specified degree will be used.
    If a vector of length 2, the two values
    indicate the degrees of trimming at the lower and upper ends of the
    \dQuote{conceptual sample} (Elamir and Seheult, 2003) of order statistics
    that is used to define the trimmed \eqn{L}-moments.
    }
}
\details{\code{samlmu} and \code{samlmu.s} are functionally identical.
\code{samlmu} calls a Fortran routine internally, and is usually faster.
\code{samlmu.s} is written entirely in the S language; it is provided
so that users can conveniently see how the calculations are done.

\code{.samlmu} is a \dQuote{bare-bones} version for use in programming.
It gives an error if \code{x} contains missing values,
computes \eqn{L}-moment ratios and not \eqn{L}-moments,
does not give a warning if all the elements of \code{x} are equal,
and returns its result in an unnamed vector.

Sample \eqn{L}-moments are defined in Hosking (1990).
Calculations use the algorithm given in Hosking (1996, p.14).

Trimmed sample \eqn{L}-moments are defined as in Hosking (2007), eq. (15)
(a small extension of Elamir and Seheult (2003), eq. (16)).
They are calculated from the untrimmed sample \eqn{L}-moments
using the recursions of Hosking (2007), eqs. (12)-(13).
}
\value{
  If \code{ratios} is \code{TRUE}, a numeric vector containing
  the \eqn{L}-moments and \eqn{L}-moment ratios,
  in the order \eqn{\ell_1}{l_1}, \eqn{\ell_2}{l_2}, \eqn{t_3}{t_3}, \eqn{t_4}{t_4}, etc.
  If \code{ratios} is \code{FALSE}, a numeric vector containing the \eqn{L}-moments
  in the order \eqn{\ell_1}{l_1}, \eqn{\ell_2}{l_2}, \eqn{l_3}{l_3}, \eqn{l_4}{l_4}, etc.
}

\note{The term \dQuote{trimmed} is used in a different sense from
its usual meaning in robust statistics.
In particular, the first trimmed \eqn{L}-moment is in general not equal to
any trimmed mean of the data sample.}

\references{
Elamir, E. A. H., and Seheult, A. H. (2003).  Trimmed L-moments.
\emph{Computational Statistics and Data Analysis}, \bold{43}, 299-314.

Hosking, J. R. M. (1990).
\eqn{L}-moments: analysis and estimation of distributions
using linear combinations of order statistics.
\emph{Journal of the Royal Statistical Society, Series B}, \bold{52}, 105-124.

Hosking, J. R. M. (1996).
Fortran routines for use with the method of \eqn{L}-moments, Version 3.
Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.

Hosking, J. R. M. (2007).  Some theory and practical uses of trimmed L-moments.
\emph{Journal of Statistical Planning and Inference}, \bold{137}, 3024-3039.
}

\author{J. R. M. Hosking \email{hosking@us.ibm.com}}

\examples{
data(airquality)
samlmu(airquality$Ozone, 6)

# Trimmed L-moment ratios
samlmu(airquality$Ozone, trim=1)
}
\keyword{univar}
