\name{lmom-package}
\docType{package}
\alias{L-moments}
\alias{Lmoments}
\alias{lmom}
\alias{lmoments}
\alias{lmom-package}
\title{The lmom package}
\description{
  \R functions for use with the method of \eqn{L}-moments
}
\details{
\eqn{L}-moments are measures of the location, scale, and shape
of probability distributions or data samples.
They are based on linear combinations of order statistics.
Hosking (1990) and Hosking and Wallis (1997, chap. 2) give expositions
of the theory of \eqn{L}-moments and \eqn{L}-moment ratios.  Hosking and Wallis
(1997, Appendix) give, for many distributions in common use, expressions
for the \eqn{L}-moments of the distributions and algorithms for estimating
the parameters of the distributions by equating sample and population
\eqn{L}-moments (the \dQuote{method of \eqn{L}-moments}).  This package contains \R functions
that should facilitate the use of \eqn{L}-moment-based methods.

For each of 13 probability distributions, the package contains
functions to evaluate the cumulative distribution function and quantile
function of the distribution, to calculate the \eqn{L}-moments given the
parameters and to calculate the parameters given the low-order
\eqn{L}-moments.  These functions are as follows.

\code{cdf...}
computes the cumulative distribution function of the distribution.

\code{qua...}
computes the quantile function (inverse cumulative distribution function)
of the distribution.

\code{lmr...}
calculates the \eqn{L}-moment ratios of the distribution given its
parameters.

\code{pel...}
calculates the parameters of the distribution given its \eqn{L}-moments.
When the \eqn{L}-moments are the sample \eqn{L}-moments of a set of data,
the resulting parameters are of course the
\dQuote{method of \eqn{L}-moments} estimates of the parameters.

Here \code{...} is a three-letter code used to identify the
distribution, as given in the table below.
For example the cumulative distribution function of the gamma distribution is
\code{cdfgam}.

\tabular{ll}{
\code{exp} \tab exponential \cr
\code{gam} \tab gamma \cr
\code{gev} \tab generalized extreme-value \cr
\code{glo} \tab generalized logistic \cr
\code{gpa} \tab generalized Pareto \cr
\code{gno} \tab generalized normal \cr
\code{gum} \tab Gumbel (extreme-value type I) \cr
\code{kap} \tab kappa \cr
\code{ln3} \tab lognormal \cr
\code{nor} \tab normal \cr
\code{pe3} \tab Pearson type III \cr
\code{wak} \tab Wakeby \cr
\code{wei} \tab Weibull \cr
}

The following functions are also contained in the package.

\code{\link{samlmu}}
  computes the sample \eqn{L}-moments of a data vector.

\code{\link{lmrp}} and \code{\link{lmrq}}
  compute the \eqn{L}-moments of a probability distribution specified
  by its cumulative distribution function (for function \code{lmrp})
  or its quantile function (for function \code{lmrq}).
  The computation uses numerical integration applied to
  a general expression for the \eqn{L}-moments of a distribution.
  Functions \code{lmrp} and \code{lmrq} can be used for any univariate
  distribution.  They are slower and usually less accurate than the
  computations carried out for specific distributions by the
  \code{lmr...} functions.

\code{\link{pelp}} and \code{\link{pelq}}
  compute the parameters of a probability distribution
  as a function of the \eqn{L}-moments.
  The computation uses function \code{lmrp} or \code{lmrq} to compute
  \eqn{L}-moments and numerical optimization to find parameter values
  for which the sample and population \eqn{L}-moments are equal.
  Functions \code{pelp} and \code{pelq} can be used for any univariate
  distribution.  They are slower and usually less accurate than the
  computations carried out for specific distributions by the
  \code{pel...} functions.

\code{\link{lmrd}}
  draws an \eqn{L}-moment ratio diagram.

\code{\link{evplot}}
  draws an \dQuote{extreme-value plot}, i.e. a quantile-quantile plot
  in which the horizontal axis is the quantile of an
  extreme-value type I (Gumbel) distribution.

\code{\link{evpoints}}, \code{\link{evdistp}}, and \code{\link{evdistq}}
  add, respectively, a set of points, a cumulative distribution function,
  and a quantile function to an extreme-value plot.
}
\section{Parameters of cumulative distribution functions and quantile functions}{
The functions \code{cdf...} (cumulative distribution functions)
and \code{qua...} (quantile functions) expect the distribution
parameters to be specified as a single vector.
This differs from the standard \R convention, in which
each parameter is a separate argument.
There are two reasons for this.
First, the single-vector parameterization is consistent with
the Fortran routines on which these \R functions are based.
Second, the single-vector parametrization is often easier to use.
For example, consider computing the 80th and 90th percentiles
of a normal distribution fitted to a set of \eqn{L}-moments
stored in a vector \code{lmom}.
In the single-vector parametrization, this is achieved by
\preformatted{
  quanor( c(.8,.9), pelnor(lmom) )
}
The separate-arguments parametrization would need a more complex
expression, such as
\preformatted{
  do.call( qnorm, c( list(.8,.9), pelnor(lmom) ) )
}

In functions (\code{lmrp}, \code{lmrq}, \code{pelp}, \code{pelq}, \code{evplot},
\code{evdistp}, \code{evdistq}) that take a cumulative distribution function
or a quantile function as an argument, the cumulative distribution function
or quantile function can use either form of parametrization.
}
\section{Relation to the LMOMENTS Fortran package}{

Functions \code{cdf...}, \code{qua...}, \code{lmr...}, \code{pel...},
and \code{samlmu} are analogous to Fortran routines from
the LMOMENTS package, version 3.04, available from StatLib at
\url{http://lib.stat.cmu.edu/general/lmoments}.
Functions \code{cdfwak} and \code{samlmu}, and all the \code{lmr...} and \code{pel...}
functions, internally call Fortran code that is derived from the
LMOMENTS package.
}
%\seealso{
%}
\references{
  Hosking, J. R. M. (1990).
  \eqn{L}-moments: analysis and estimation of distributions
  using linear combinations of order statistics.
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{52}, 105-124.

  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on L-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
\keyword{package}
