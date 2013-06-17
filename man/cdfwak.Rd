\name{cdfwak}
\alias{cdfwak}
\alias{quawak}
\title{Wakeby distribution}
\description{Distribution function and quantile function
of the Wakeby distribution.}
\usage{
cdfwak(x, para = c(0, 1, 0, 0, 0))
quawak(f, para = c(0, 1, 0, 0, 0))
}

\arguments{
\item{x}{Vector of quantiles.}
\item{f}{Vector of probabilities.}
\item{para}{Numeric vector containing the parameters of the distribution,
in the order
\eqn{\xi, \alpha, \beta, \gamma, \delta}{xi, alpha, beta, gamma, delta}.}
}

\details{The Wakeby distribution with
parameters \eqn{\xi}{xi},
\eqn{\alpha}{alpha},
\eqn{\beta}{beta},
\eqn{\gamma}{gamma} and
\eqn{\delta}{delta}
has quantile function
\deqn{x(F)=\xi+{\alpha\over\beta}\lbrace1-(1-F)^\beta\rbrace-{\gamma\over\delta}\lbrace1-(1-F)^{-\delta}\rbrace.}{x(F) = xi + alpha {1-(1-F)^beta}/beta - gamma {1-(1-F)^(-delta)}/delta .}

The parameters are restricted as in Hosking and Wallis (1997, Appendix A.11):
\itemize{
\item either \eqn{\beta+\delta>0}{beta + delta > 0} or
\eqn{\beta=\gamma=\delta=0}{beta = gamma = delta = 0};
\item if \eqn{\alpha=0}{alpha = 0} then \eqn{\beta=0}{beta = 0};
\item if \eqn{\gamma=0}{gamma = 0} then \eqn{\delta=0}{delta = 0};
\item \eqn{\gamma\ge0}{gamma >= 0};
\item \eqn{\alpha+\gamma\ge0}{alpha + gamma >= 0}.
}

The distribution has a lower bound at \eqn{\xi}{xi} and,
if \eqn{\delta<0}{delta<0}, an upper bound at
\eqn{\xi+\alpha/\beta-\gamma/\delta}{xi+alpha/beta-gamma/delta}.


The generalized Pareto distribution is the special case
\eqn{\alpha=0}{alpha=0} or \eqn{\gamma=0}{gamma=0}.
The exponential distribution is the special case
\eqn{\beta=\gamma=\delta=0}{beta=gamma=delta=0}.
The uniform distribution is the special case
\eqn{\beta=1}{beta=1}, \eqn{\gamma=\delta=0}{gamma=delta=0}.
}

\value{\code{cdfwak} gives the distribution function;
       \code{quawak} gives the quantile function.}


\note{The functions expect the distribution parameters in a vector,
rather than as separate arguments as in the standard \R
distribution functions \code{pnorm}, \code{qnorm}, etc.}


\seealso{
\code{\link{cdfgpa}} for the generalized Pareto distribution.

\code{\link{cdfexp}} for the exponential distribution.
}

\references{
  Hosking, J. R. M. and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on L-moments},
  Cambridge University Press, Appendix A.11.
}

\examples{
# Random sample from the Wakeby distribution
# with parameters xi=0, alpha=30, beta=20, gamma=1, delta=0.3.
quawak(runif(100), c(0,30,20,1,0.3))
}
\keyword{distribution}
