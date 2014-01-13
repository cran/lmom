\name{lmrp}
\alias{lmrp}
\alias{lmrq}
\title{L-moments of a general probability distribution}
\description{Computes the \eqn{L}-moments or trimmed \eqn{L}-moments
of a probability distribution
given its cumulative distribution function (for function \code{lmrp})
or quantile function (for function \code{lmrq}).
}
\usage{
lmrp(pfunc, ..., bounds=c(-Inf,Inf), symm=FALSE, order=1:4,
     ratios=TRUE, trim=0, acc=1e-6, subdiv=100, verbose=FALSE)

lmrq(qfunc, ..., symm=FALSE, order=1:4, ratios=TRUE, trim=0,
     acc=1e-6, subdiv=100, verbose=FALSE)
}

\arguments{
  \item{pfunc}{Cumulative distribution function.}
  \item{qfunc}{Quantile function.}
  \item{\dots}{Arguments to \code{pfunc} or \code{qfunc}.}
  \item{bounds}{Either a vector of length 2, containing the lower
    and upper bounds of the distribution, or a function that calculates
    these bounds given the distribution parameters as inputs.}
  \item{symm}{
    For \code{lmrq}, a logical value indicating whether
    the distribution is symmetric about its median.

    For \code{lmrp}, either the logical value \code{FALSE}
    to indicate that the distribution is not symmetric,
    or a numeric value to indicate that the distribution is symmetric
    and that the specified value is the center of symmetry.

    If the distribution is symmetric, odd-order \eqn{L}-moments
    are exactly zero and the symmetry is used to slightly speed up
    the computation of even-order \eqn{L}-moments.}
  \item{order}{Orders of the \eqn{L}-moments and \eqn{L}-moment ratios to be computed.}
  \item{ratios}{Logical. If \code{FALSE}, \eqn{L}-moments are computed;
    if \code{TRUE} (the default), \eqn{L}-moment ratios are computed.}
  \item{trim}{Degree of trimming.
    If a single value, symmetric trimming of the specified degree will be used.
    If a vector of length 2, the two values
    indicate the degrees of trimming at the lower and upper ends of the
    \dQuote{conceptual sample} (Elamir and Seheult, 2003) of order statistics
    that is used to define the trimmed \eqn{L}-moments.}
  \item{acc}{Requested accuracy.  The function will try to achieve
    this level of accuracy, as relative error for \eqn{L}-moments and absolute
    error for \eqn{L}-moment ratios.}
  \item{subdiv}{Maximum number of subintervals used in numerical integration.}
  \item{verbose}{Logical.  If \code{FALSE}, only the values of the \eqn{L}-moments
    and \eqn{L}-moment ratios are returned.  If \code{TRUE}, more details of the
    numerical integration are returned: see \dQuote{Value} section below.}
}

\details{Computations use expressions in Hosking (2007):
eq. (7) for \code{lmrp}, eq. (5) for \code{lmrq}.
Integrals in those expressions are computed by numerical integration.
}

\value{If \code{verbose=FALSE} and \code{ratios=FALSE},
a numeric vector containing the \eqn{L}-moments.

If \code{verbose=FALSE} and \code{ratios=TRUE},
a numeric vector containing the \eqn{L}-moments (of orders 1 and 2)
and \eqn{L}-moment ratios (of orders 3 and higher).

If \code{verbose=TRUE}, a data frame with columns as follows:
\item{value}{\eqn{L}-moments (if \code{ratios=FALSE}),
or \eqn{L}-moments and \eqn{L}-moment ratios (if \code{ratios=TRUE}).}
\item{abs.error}{Estimate of the absolute error in the computed value.}
\item{message}{\code{"OK"} or a character string giving the error message
resulting from the numerical integration.}
}

\section{Arguments of cumulative distribution functions and quantile functions}{
\code{pfunc} and \code{qfunc} can be either the standard \R form of
cumulative distribution function or quantile function
(i.e. for a distribution with \eqn{r} parameters, the first argument is the
variate \eqn{x} or the probability \eqn{p} and the next \eqn{r} arguments
are the parameters of the distribution) or the \code{cdf\dots} or
\code{qua\dots} forms used throughout the \pkg{lmom} package
(i.e. the first argument is the variate \eqn{x} or probability \eqn{p}
and the second argument is a vector containing the parameter values).
Even for the \R form, however, starting values for the parameters
are supplied as a vector \code{start}.

If \code{bounds} is a function, its arguments must match
the distribution parameter arguments of \code{pfunc}:
either a single vector, or a separate argument for each parameter.
}

\note{
In package \pkg{lmom} versions 1.6 and earlier, the \dQuote{Details} section stated that
\dQuote{Integrals in those expressions are computed by numerical integration,
using the \R function \code{integrate}}.
As of version 2.0, numerical integration uses an internal function
that directly calls (slightly modified versions of)
the Fortran routines in QUADPACK (Piessens et al. 1983).
\R{}'s own \code{integrate} function uses C code \dQuote{based on} the
QUADPACK routines, but in \R{} versions 2.12.0 through 3.0.1
did not in every case reproduce the results that
would have been obtained with the Fortran code
(this is \R{} bug PR#15219).
}

\references{
  Elamir, E. A. H., and Seheult, A. H. (2003).  Trimmed L-moments.
  \emph{Computational Statistics and Data Analysis}, \bold{43}, 299-314.

  Hosking, J. R. M. (2007).  Some theory and practical uses of trimmed L-moments.
  \emph{Journal of Statistical Planning and Inference}, \bold{137}, 3024-3039.

  Piessens, R., deDoncker-Kapenga, E., Uberhuber, C., and Kahaner, D. (1983).
  \emph{Quadpack: a Subroutine Package for Automatic Integration}.
  Springer Verlag.
}
\author{J. R. M. Hosking \email{hosking@us.ibm.com}}

\section{Warning}{Arguments \code{bounds}, \code{symm}, \code{order},
\code{ratios}, \code{trim}, \code{acc}, \code{subdiv}, and \code{verbose}
cannot be abbreviated and must be specified by their full names
(if abbreviated, the names would be matched to the arguments of
\code{pfunc} or \code{qfunc}).}

\seealso{
\code{\link{lmrexp}} to compute (untrimmed) \eqn{L}-moments of specific distributions.

\code{\link{samlmu}} to compute (trimmed or untrimmed) \eqn{L}-moments of a data sample.

\code{\link{pelp}} and \code{\link{pelexp}},
to compute the parameters of a distribution given its (trimmed or untrimmed) \eqn{L}-moments.
}

\examples{
## Generalized extreme-value (GEV) distribution
## - three ways to get its L-moments
lmrp(cdfgev, c(2,3,-0.2))
lmrq(quagev, c(2,3,-0.2))
lmrgev(c(2,3,-0.2), nmom=4)

## GEV bounds specified as a vector
##
lmrp(cdfgev, c(2,3,-0.2), bounds=c(-13,Inf))
##
## GEV bounds specified as a function -- single vector of parameters
##
gevbounds <- function(para) {
  k <- para[3]
  b <- para[1]+para[2]/k
  c(ifelse(k<0, b,- Inf), ifelse(k>0, b, Inf))
}
lmrp(cdfgev, c(2,3,-0.2), bounds=gevbounds)
##
## GEV bounds specified as a function -- separate parameters
##
pgev <- function(x, xi, alpha ,k)
  pmin(1, pmax(0, exp(-((1-k*(x-xi)/alpha)^(1/k)))))
pgevbounds <- function(xi,alpha,k) {
  b <- xi+alpha/k
  c(ifelse(k<0, b, -Inf), ifelse(k>0, b, Inf))
}
lmrp(pgev, xi=2, alpha=3, k=-0.2, bounds=pgevbounds)


## Normal distribution
lmrp(pnorm)
lmrp(pnorm, symm=0)
lmrp(pnorm, mean=2, sd=3, symm=2)
# For comparison, the exact values
lmrnor(c(2,3), nmom=4)

# Many L-moment ratios of the exponential distribution
# This may warn that "the integral is probably divergent"
lmrq(qexp, order=3:20)

# ... nonetheless the computed values seem accurate:
# compare with the exact values, tau_r = 2/(r*(r-1)):
cbind(exact=2/(3:20)/(2:19), lmrq(qexp, order=3:20, verbose=TRUE))

# Of course, sometimes the integral really is divergent
\dontrun{
lmrq(function(p) (1-p)^(-1.5))
}

# And sometimes the integral is divergent but that's not what
# the warning says (at least on the author's system)
lmrp(pcauchy)

# Trimmed L-moments for Cauchy distribution are finite
lmrp(pcauchy, symm=0, trim=1)

# Works for discrete distributions too, but often requires
# a larger-than-default value of 'subdiv'
lmrp(ppois, lambda=5, subdiv=1000)
}
\keyword{distribution}
