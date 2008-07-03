\name{cdfpe3}
\alias{cdfpe3}
\alias{quape3}
\title{Pearson type III distribution}
\description{Distribution function and quantile function
of the Pearson type III distribution}
\usage{
cdfpe3(x, para = c(0, 1, 0))
quape3(f, para = c(0, 1, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\mu, \sigma, \gamma}{mu, sigma, gamma}
(location, scale, shape).}
}

\details{The Pearson type III distribution contains as special cases
the usual three-parameter gamma distribution
(a shifted version of the gamma distribution)
with a finite lower bound and positive skewness;
the normal distribution,
and the reverse three-parameter gamma distribution,
with a finite upper bound and negative skewness.
The distribution's parameters are the first three (ordinary) moment ratios:
\eqn{\mu}{mu} (the mean, a location parameter),
\eqn{\sigma}{sigma} (the standard deviation, a scale parameter) and
\eqn{\gamma}{gamma} (the skewness, a shape parameter).

If \eqn{\gamma\ne0}{gamma != 0}, let \eqn{\alpha=4/\gamma^2}{alpha=4/gamma^2},
\eqn{\beta={\scriptstyle 1 \over \scriptstyle 2}\sigma|\gamma|}{(1/2) sigma |gamma|},
\eqn{\xi=\mu-2\sigma/\gamma}{xi = mu - 2 sigma / gamma}.
The probability density function is
\deqn{f(x)={|x-\xi|^{\alpha-1}\exp(-|x-\xi|/\beta) \over \beta^\alpha\Gamma(\alpha)}}{f(x) = |x-xi|^(\alpha-1) exp(-|x-\xi|/\beta) / ( beta^alpha Gamma(alpha) )}
with \eqn{x}{x} bounded by \eqn{\xi}{xi} from below if \eqn{\gamma>0}{gamma>0}
and from above if \eqn{\gamma<0}{gamma<0}.
If \eqn{\gamma=0}{gamma=0}, the distribution ius a normal distribution
with mean \eqn{\mu}{mu} and standard deviation \eqn{\sigma}{sigma}.

The Pearson type III distribution is usually regarded as consisting of
just the case \eqn{\gamma>0}{gamma>0} given above, and is usually
parametrized by \eqn{\alpha}{alpha}, \eqn{\beta}{beta} and \eqn{\xi}{xi}.
Our parametrization extends the distribution to include
the usual Pearson type III distributions,
with positive skewness and lower bound \eqn{\xi}{xi},
reverse Pearson type III distributions,
with negative skewness and upper bound \eqn{\xi}{xi},
and the Normal distribution, which is included as a special
case of the distribution rather than as the unattainable limit
\eqn{\alpha\rightarrow\infty}{alpha -> infinity}.
This enables the Pearson type III distribution to be used when the skewness of
the observed data may be negative.
The parameters \eqn{\mu}{mu}, \eqn{\sigma}{sigma} and \eqn{\gamma}{gamma}
are the conventional moments of the distribution.

The gamma distribution is obtained when \eqn{\gamma>0}{gamma>0}
and \eqn{\mu=2\sigma/\gamma}{mu = 2 sigma / gamma }.
The normal distribution is the special case \eqn{\gamma=0}{gamma=0}.
The exponential distribution is the special case \eqn{\gamma=2}{gamma=2}.
}

\value{\code{cdfpe3} gives the distribution function;
       \code{quape3} gives the quantile function.}

\references{
  Hosking, J. R. M. and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on L-moments},
  Cambridge University Press, Appendix A.10.
}

\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgam}} for the gamma distribution.

\code{\link{cdfnor}} for the normal distribution.
}

\examples{
# Random sample from the Pearson type III distribution
# with parameters mu=1, alpha=2, gamma=3.
quape3(runif(100), c(1,2,3))

# The Pearson type III distribution with parameters
# mu=12, sigma=6, gamma=1, is the gamma distribution
# with parameters alpha=4, beta=3.  An illustration:
fval<-seq(0.1,0.9,by=0.1)
cbind(fval, qgamma(fval, shape=4, scale=3), quape3(fval, c(12,6,1)))
}
\keyword{distribution}
