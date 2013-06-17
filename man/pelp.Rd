\name{pelp}
\alias{pelp}
\alias{pelq}
\title{Parameter estimation for a general distribution by the method of L-moments}
\description{Computes the parameters of a probability distribution
as a function of the \eqn{L}-moments or trimmed \eqn{L}-moments.
}
\usage{
pelp(lmom, pfunc, start, bounds = c(-Inf, Inf),
     type = c("n", "s", "ls", "lss"),
     ratios = NULL, trim = NULL, method = "nlm", acc = 1e-5, ...)

pelq(lmom, qfunc, start, type = c("n", "s", "ls", "lss"),
     ratios = NULL, trim = NULL, method = "nlm", acc = 1e-5, ...)
}

\arguments{
  \item{lmom}{Numeric vector containing the \eqn{L}-moments of the distribution
    or of a data sample.}
  \item{pfunc}{Cumulative distribution function of the distribution.}
  \item{qfunc}{Quantile function of the distribution.}
  \item{start}{Vector of starting values for the parameters.}
  \item{bounds}{Either a vector of length 2, containing the lower
    and upper bounds of the distribution, or a function that calculates
    these bounds given the distribution parameters as inputs.}
  \item{type}{Type of distribution, i.e. how it is parametrized.
    Must be one of the following:
    \describe{
    \item{\code{"ls"}}{The distribution has a location parameter and a scale parameter.}
    \item{\code{"lss"}}{The distribution has a location parameter and a scale parameter,
      and is symmetric about its median.}
    \item{\code{"s"}}{The distribution has a scale parameter but not a location parameter.}
    \item{\code{"n"}}{The distribution has neither a location parameter
      nor a scale parameter.}}

    For more details, see the \dQuote{Distribution type} section below.
  }
  \item{ratios}{Logical or \code{NULL}.
    If \code{FALSE}, \code{lmom} should contain \eqn{L}-moments;
    if \code{TRUE}, \code{lmom} should contain \eqn{L}-moment ratios.
    If \code{NULL} and \code{lmom} has names, the contents of \code{lmom}
    will be inferred from these names - see section
    \dQuote{Inferring \sQuote{ratios} and \sQuote{trim}} below.
    The default value (if \code{ratios} is \code{NULL} and \code{lmom} has no names)
    is \code{TRUE}.}
  \item{trim}{The degree of trimming corresponding to the \eqn{L}-moments in \code{lmom}.
    Can be a single value or a vector length 2, as for \code{samlmu}.
    Can also be \code{NULL}: in this case if \code{lmom} has names,
    the degree of trimming will be inferred from these names - see section
    \dQuote{Inferring \sQuote{ratios} and \sQuote{trim}} below.
    The default value (if \code{trim} is \code{NULL} and \code{lmom} has no names)
    is \code{0}.}
  \item{method}{Method used to estimate the shape parameters
    (i.e. all parameters other than the location and scale parameters, if any).
    Valid values are \code{"nlm"} (the default), \code{"uniroot"}
    (which is valid only if the distribution has at most one shape parameter),
    and any of the values of the \code{method} argument of function \code{optim}.
    See the \dQuote{Details} section below.
  }
  \item{acc}{Requested accuracy for the estimated parameters.
    This will be absolute accuracy for shape parameters,
    relative accuracy for a scale parameter, and
    absolute accuracy of the location parameter divided by the scale parameter
    for a location parameter.}
  \item{\dots}{Further arguments will be passed to the optimization function
    (\code{nlm}, \code{uniroot}, or \code{optim}).}
}

\details{For shape parameters, numerical optimization is used to
find parameter values for which the population \eqn{L}-moments or \eqn{L}-moment ratios
are equal to the values supplied in \code{lmom}.
Computation of \eqn{L}-moments or \eqn{L}-moment ratios
uses functions \code{lmrp} (for \code{pelp}) or \code{lmrq} (for \code{pelq}).
Numerical optimization uses \R functions \code{nlm} (if \code{method="nlm"}),
\code{uniroot} (if \code{method="uniroot"}), or
\code{optim} with the specified method (for the other values of \code{method}).
Function \code{uniroot} uses one-dimensional root-finding,
while functions \code{nlm} and \code{optim} try to minimize
a criterion function that is the sum of squared differences between the
population \eqn{L}-moments or \eqn{L}-moment ratios and the values supplied in \code{lmom}.
Location and scale parameters are then estimated noniteratively.
In all cases, the calculation of population \eqn{L}-moments and
\eqn{L}-moment ratios is performed by function \code{lmrp} or \code{lmrq}
(when using \code{pelp} or \code{pelq} respectively).

This approach is very crude.  Nonetheless, it is often effective in practice.
As in all numerical optimizations, success may depend on the way that
the distribution is parametrized and on the particular choice of
starting values for the parameters.
}

\value{A list with components:
\item{para}{Numeric vector containing the estimated parameters
of the distribution.}
\item{code}{An integer indicating the result of the numerical optimization
used to estimate the shape parameters.  It is \code{0} if there
are no shape parameters.  In general, values \code{1} and \code{2}
indicate successful convergence of the iterative procedure,
a value of \code{3} indicates that the iteration may not have converged,
and values of \code{4} or more indicate that the iteration did not converge.
Specifically, \code{code} is:

For method \code{"nlm"}, the \code{code} component
of the return value from \code{nlm}.

For method \code{"uniroot"}, \code{1} if the estimated
precision of the shape parameter is less than or equal to \code{acc},
and \code{4} otherwise.

For the other methods, the \code{convergence} component
of the return value from \code{optim}.
}
}

\section{Further details of arguments}{
The length of \code{lmom} should be (at least) the highest
order of \eqn{L}-moment used in the estimation procedure.  For a distribution
with \eqn{r} parameters this is
\eqn{2r-2} if \code{type="lss"} and \eqn{r} otherwise.

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

It is assumed that location and scale parameters come first in
the set of parameters of the distribution. Specifically:
if \code{type="ls"} or \code{type="lss"}, it is assumed
that the first parameter is the location parameter and
that the second parameter is the scale parameter;
if \code{type="s"} it is assumed
that the first parameter is the scale parameter.

It is important that the length of \code{start} be equal to the number
of parameters of the distribution.  Starting values for
location and scale parameters should be included in \code{start},
even though they are not used.
If \code{start} has the wrong length, it is possible that
meaningless results will be returned without any warning being issued.

}

\section{Distribution type }{
The \code{type} argument affects estimation as follows.
We assume that location and scale parameters are \eqn{\xi} and \eqn{\alpha}
respectively, and that the shape parameters (if there are any)
are collectively designated by \eqn{\theta}.

If \code{type="ls"}, then the \eqn{L}-moment ratios \eqn{\tau_3, \tau_4, \ldots}
depend only on the shape parameters. If there are any shape parameters,
they are estimated by equating the sample \eqn{L}-moment ratios of orders
3, 4, etc., to the population \eqn{L}-moment ratios
and solving the resulting equations for the shape parameters
(using as many equations as there are shape parameters).
The \eqn{L}-moment \eqn{\lambda_2} is a multiple of \eqn{\alpha}, the multiplier
being a function only of \eqn{\theta}.
\eqn{\alpha} is estimated by dividing the second sample \eqn{L}-moment
by the multiplier function evaluated at the estimated value of \eqn{\theta}.
The \eqn{L}-moment \eqn{\lambda_1} is \eqn{\xi} plus
a function of \eqn{\alpha} and \eqn{\theta}.
\eqn{\xi} is estimated by subtracting from the first sample \eqn{L}-moment
the function evaluated at the estimated values of
\eqn{\alpha} and \eqn{\theta}.

If \code{type="lss"}, then
the \eqn{L}-moment ratios of odd order, \eqn{\tau_3, \tau_5, \ldots}, are zero and
the \eqn{L}-moment ratios of even order, \eqn{\tau_4, \tau_6, \ldots},
depend only on the shape parameters. If there are any shape parameters,
they are estimated by equating the sample \eqn{L}-moment ratios of orders
4, 6, etc., to the population \eqn{L}-moment ratios
and solving the resulting equations for the shape parameters
(using as many equations as there are shape parameters).
Parameters \eqn{\alpha} and \eqn{\xi} are estimated as in case when \code{type="ls"}.

If \code{type="s"}, then the \eqn{L}-moments divided by \eqn{\lambda_1},
i.e. \eqn{\lambda_2/\lambda_1, \lambda_3/\lambda_1, \ldots},
depend only on the shape parameters. If there are any shape parameters,
they are estimated by equating the sample \eqn{L}-moments
(divided by the first sample \eqn{L}-moment) of orders 2, 3, etc.,
to the corresponding population \eqn{L}-moments
(divided by the first population \eqn{L}-moment)
and solving the resulting equations
(as many equations as there are shape parameters).
The \eqn{L}-moment \eqn{\lambda_1} is a multiple of \eqn{\alpha}, the multiplier
being a function only of \eqn{\theta}.
\eqn{\alpha} is estimated by dividing the first sample \eqn{L}-moment
by the multiplier function evaluated at the estimated value of \eqn{\theta}.

If \code{type="n"}, then all parameters are shape parameters.
They are estimated by equating the sample \eqn{L}-moments of orders
1, 2, etc., to the population \eqn{L}-moments
and solving the resulting equations for the parameters
(using as many equations as there are parameters).

}

\section{Inferring \sQuote{ratios} and \sQuote{trim}}{
If \code{ratios} or \code{trim} is \code{NULL}, appropriate values will be inferred
by inspecting the names of \code{lmom}. It is assumed that \code{lmom}
was generated by a call to \code{samlmu}, \code{lmrp}, or \code{lmrq};
in this case its names will reflect the values of \code{ratios} and \code{trim}
used in that  call.
If in this case \code{lmom} has no names, default values
\code{ratios=TRUE} and \code{trim=0} will be used.

This inference is made in order to reduce the need to specify the
orders of trimming repetitivelty.
For example, a distribution with quantile function \code{qfunc} can be
fitted to (1,1)-trimmed \eqn{L}-moments of data in \code{x} by
\preformatted{
  lmom <- samlmu(x, trim=1)
  fit <- pelq(lmom, qfunc, start=...)
}
There is no need to specify \code{trim} both in the call to \code{samlmu}
and the call to \code{pelq}.
}

\author{J. R. M. Hosking \email{hosking@us.ibm.com}}

\seealso{
\code{\link{pelexp}} for parameter estimation of specific distributions.
}

\examples{
## Gamma distribution -- rewritten so that its first parameter
## is a scale parameter
my.pgamma <- function(x, scale, shape) pgamma(x, shape=shape, scale=scale)
pelp(c(5,2), my.pgamma, start=c(1,1), bounds=c(0,Inf), type="s")
# We can also do the estimation suppressing our knowledge
# that one parameter is a shape parameter.
pelp(c(5,2), my.pgamma, start=c(1,1), bounds=c(0,Inf), type="n")
rm(my.pgamma)

## Kappa distribution -- has location, scale and 2 shape parameters
# Estimate via pelq
pel.out <- pelq(c(10,5,0.3,0.15), quakap, start=c(0,1,0,0), type="ls")
pel.out
# Check that L-moments of estimated distribution agree with the
# L-moments input to pelq()
lmrkap(pel.out$para)
# Compare with the distribution-specific routine pelkap
pelkap(c(10,5,0.3,0.15))
rm(pel.out)

# Similar results -- what's the advantage of the specific routine?
system.time(pelq(c(10,5,0.3,0.15), quakap, start=c(0,1,0,0), type="ls"))
system.time(pelkap(c(10,5,0.3,0.15)))

# Caution -- pelq() will not check that estimates are reasonable
lmom <- c(10,5,0.2,0.25)
pel.out <- pelq(lmom, quakap, start=c(0,1,0,0), type="ls")
pel.out
lmrkap(pel.out$para) # should be close to lmom, but tau_3 and tau_4 are not
# What happened? pelkap will tell us
try(pelkap(lmom))
rm(lmom, pel.out)

## Inverse Gaussian -- don't have explicit estimators for this
## distribution, but can use numerical methods
#
# CDF of inverse gaussian distribution
pig <- function(x, mu, lambda) {
  temp <- suppressWarnings(sqrt(lambda/x))
  xx <- pnorm(temp*(x/mu-1))+exp(2*lambda/mu+pnorm(temp*(x/mu+1),
              lower.tail=FALSE, log.p=TRUE))
  out <- ifelse(x<=0, 0, xx)
  out
}
# Fit to ozone data
data(airquality)
(lmom<-samlmu(airquality$Ozone))
pel.out <- pelp(lmom[1:2], pig, start=c(10,10), bounds=c(0,Inf))
pel.out
# First four L-moments of fitted distribution,
# for comparison with sample L-moments
lmrp(pig, pel.out$para[1], pel.out$para[2], bounds=c(0,Inf))
rm(pel.out)

## A Student t distribution with location and scale parameters
#
qstu <- function(p, xi, alpha, df) xi + alpha * qt(p, df)
# Estimate parameters.  Distribution is symmetric: use type="lss"
pelq(c(3,5,0,0.2345), qstu, start=c(0,1,10), type="lss")
# Doesn't converge (at least on the author's system) --
# try a different parametrization
qstu2 <- function(p, xi, alpha, shape) xi + alpha * qt(p, 1/shape)
# Now it converges
pelq(c(3,5,0,0.2345), qstu2, start=c(0,1,0.1), type="lss")
# Or try a different optimization method
pelq(c(3,5,0,0.2345), qstu, start=c(0,1,10), type="lss",
    method="uniroot", lower=2, upper=100)

## With trimmed L-moments, we can fit this distribution even when
## it does not have a finite mean ('df' less than 1)
set.seed(123456)
dat <- qstu(runif(1000), xi=3, alpha=5, df=0.75)
lmom <- samlmu(dat, trim=1)
lmom
# Note that pelq() infers 'trim=1' from the names of 'lmom'
pelq(lmom, qstu, start=c(0,1,10), type="lss",  method="uniroot",
  lower=0.51, upper=100)

rm(qstu, qstu2, dat, lmom)

}
\keyword{models}
