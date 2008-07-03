\name{cdfkap}
\alias{cdfkap}
\alias{quakap}
\title{Kappa distribution}
\description{Distribution function and quantile function
of the kappa distribution.}
\usage{
cdfkap(x, para = c(0, 1, 0, 0))
quakap(f, para = c(0, 1, 0, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order \eqn{\xi, \alpha, k, h}{xi, alpha, k, h} (location, scale,
shape, shape).}
}

\details{The kappa distribution with
location parameter \eqn{\xi}{xi},
scale parameter \eqn{\alpha}{alpha} and
shape parameters \eqn{k}{k} and \eqn{h}{h}
has quantile function
\deqn{x(F)=\xi-{\alpha\over k}(1-({1-F^h \over h})^k).}{x(F) = xi - alpha (1 - ( (1-F^h)/h )^k) / k .}

Its special cases include the
generalized logistic (\eqn{h=-1}{h=-1}),
generalized extreme-value (\eqn{h=0}{h=0}),
generalized Pareto (\eqn{h=1}{h=1}),
logistic (\eqn{k=0}{k=0}, \eqn{h=-1}{h=-1}),
Gumbel (\eqn{k=0}{k=0}, \eqn{h=0}{h=0}),
exponential (\eqn{k=0}{k=0}, \eqn{h=1}{h=1}), and
uniform (\eqn{k=1}{k=1}, \eqn{h=1}{h=1}) distributions.
}

\value{\code{cdfkap} gives the distribution function;
       \code{quakap} gives the quantile function.}

\references{
  Hosking, J. R. M. (1994). The four-parameter kappa distribution.
  \emph{IBM Journal of Research and Development}, \bold{38}, 251-258.

  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on L-moments},
  Cambridge University Press, Appendix A.10.
}

\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfglo}} for the generalized logistic distribution,
\code{\link{cdfgev}} for the generalized extreme-value distribution,
\code{\link{cdfgpa}} for the generalized Pareto distribution,
\code{\link{cdfgum}} for the Gumbel distribution,
\code{\link{cdfexp}} for the exponential distribution.
}

\examples{
# Random sample from the kappa distribution
# with parameters xi=0, alpha=1, k=-0.5, h=0.25.
quakap(runif(100), c(0,1,-0.5,0.25))
}
\keyword{distribution}
