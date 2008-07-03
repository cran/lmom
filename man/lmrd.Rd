\name{lmrd}
\alias{lmrd}
\title{L-moment ratio diagram}
\description{Draws an \eqn{L}-moment ratio diagram.}
\usage{
lmrd(x, y, distributions = "GLO GEV GPA GNO PE3", twopar,
     xlim, ylim, pch=3, col, lty, lwd=1,
     legend.lmrd = TRUE, xlegend, ylegend,
     xlab = expression(italic(L) * "-skewness"),
     ylab = expression(italic(L) * "-kurtosis"), ...)
}

\arguments{
\item{x}{Numeric vector of \eqn{L}-skewness values.

Alternatively, if the \code{y} argument is omitted, \code{x} can be a vector
of \eqn{L}-moment ratios, as returned for example by function \code{samlmu}.
In this case the third and fourth elements of \code{x} are taken to be
\eqn{L}-skewness and \eqn{L}-kurtosis respectively and are plotted.}

\item{y}{Numeric vector of \eqn{L}-kurtosis values.}

\item{distributions}{Indicates the three-parameter distributions
whose \eqn{L}-skewness--\eqn{L}-kurtosis relations
are to be plotted as lines on the diagram.
The following distributions are recognized:
\tabular{ll}{
\code{GLO} \tab generalized logistic \cr
\code{GEV} \tab generalized extreme-value \cr
\code{GPA} \tab generalized Pareto \cr
\code{GNO} \tab generalized normal (lognormal) \cr
\code{PE3} \tab Pearson type III \cr
\code{WAK.LB} \tab lower bound of \eqn{L}-kurtosis for given \eqn{L}-skewness, \cr
              \tab for the Wakeby distribution. \cr
\code{ALL.LB} \tab lower bound of \eqn{L}-kurtosis for given \eqn{L}-skewness, \cr
              \tab for all distributions.
}

% Alternative to a table:
%\itemize{
%\item GLO -- generalized logistic
%\item GEV -- generalized extreme-value
%\item GPA -- generalized Pareto
%\item GNO -- generalized normal (lognormal)
%\item PE3 -- Pearson type III
%\item WAK.LB -- lower bound of \eqn{L}-kurtosis for given \eqn{L}-skewness, for the Wakeby distribution.
%\item ALL.LB -- lower bound of \eqn{L}-kurtosis for given \eqn{L}-skewness, for all distributions.
%}

The argument should be either a character vector each of whose elements
is one of the above abbreviations or a character string
containing one or more of the abbreviations separated by blanks.
The specified \eqn{L}-skewness--\eqn{L}-kurtosis curves will be plotted.

If no three-parameter distributions are to be plotted,
specify \code{distributions} to be \code{FALSE} or
the empty string, \code{""}.
}
\item{twopar}{
Two-parameter distributions whose (\eqn{L}-skewness, \eqn{L}-kurtosis) values
are to be plotted as points on the diagram.
The following distributions are recognized:
\tabular{ll}{
\code{E} or \code{EXP} \tab exponential \cr
\code{G} or \code{GUM} \tab Gumbel \cr
\code{L} or \code{LOG} \tab logistic \cr
\code{N} or \code{NOR} \tab normal \cr
\code{U} or \code{UNI} \tab uniform
}

The argument should be either a character vector each of whose elements
is one of the above abbreviations or a character string
containing one or more of the abbreviations separated by blanks.
\eqn{L}-skewness--\eqn{L}-kurtosis points for the specified distributions will be plotted
and given one-character labels.

The default is to plot the two-parameter distributions that are special
cases of the three-parameter distributions specified in the
\code{distributions} parameter.  Thus for example if
\code{distributions="GPA PE3"}, the default for \code{twopar}
is \code{"EXP NOR UNI"}: NOR is a special case of PE3, UNI of GPA,
EXP of both GPA and PE3.

If no two-parameter distributions are to be plotted,
specify \code{twopar} to be \code{FALSE} or the empty string, \code{""}.
}

\item{xlim}{x axis limits.
Default: \code{c(0, 0.6)}, expanded if necessary to cover the range of the data.}

\item{ylim}{y axis limits.
Default: \code{c(0, 0.4)}, expanded if necessary to cover the range of the data.}

\item{pch}{Plotting character to be used for the plotted
(\eqn{L}-skewness, \eqn{L}-kurtosis) points.}

\item{col}{Vector specifying the colors to be used for the lines on the plot.
For the default values, see the description of the \code{lty} argument below.}

\item{lty}{Vector specifying the line types to be used for the lines on the plot.

By default, colors and line types are matched to the distributions given
in the \code{distributions} parameter, as follows:
\tabular{ll}{
% S-plus colors
%GLO \tab color  6, line type 1 (blue, solid) \cr
%GEV \tab color  4, line type 1 (green, solid) \cr
%GPA \tab color 15, line type 1 (red, solid) \cr
%GNO \tab color  1, line type 1 (black, solid) \cr
%PE3 \tab color 13, line type 1 (cyan, solid) \cr
%WAK.LB \tab color 15, line type 4 (red, dashed) \cr
%ALL.LB \tab color  1, line type 4 (black, dashed)
% R colors
\code{GLO}    \tab color \code{"blue"},  solid line \cr
\code{GEV}    \tab color \code{"green"}, solid line \cr
\code{GPA}    \tab color \code{"red"},   solid line \cr
\code{GNO}    \tab color \code{"black"}, solid line \cr
\code{PE3}    \tab color \code{"cyan"},  solid line \cr
\code{WAK.LB} \tab color \code{"red"},   dashed line\cr
\code{ALL.LB} \tab color \code{"black"}, dashed line
}
}

\item{lwd}{Vector specifying the line widths to be used for the lines on the plot.}

\item{legend.lmrd}{Controls whether a legend,
identifying the \eqn{L}-skewness--\eqn{L}-kurtosis
relations of the three-parameter distributions, is plotted.
Not used if \code{distributions=FALSE}.}

\item{xlegend}{x coordinate of the upper left corner of the legend.
Default: the minimum x value.
Not used if \code{distributions=FALSE} or \code{legend.lmrd=FALSE}.}

\item{ylegend}{y coordinate of the upper left corner of the legend.
Default: the maximum y value.
Not used if \code{distributions=FALSE} or \code{legend.lmrd=FALSE}.}

\item{xlab}{X axis label.}

\item{ylab}{Y axis label.}

\item{\dots}{Additional arguments are passed to the \code{matplot} function
that draws the axis box and the lines for three-parameter distributions.}
}
\details{
\code{lmrd} calls a sequence of graphics functions:
\code{matplot} for the axis box and the curves for three-parameter distributions;
\code{points} for the points for two-parameter distributions and
\code{text} for their labels; \code{legend} for the legend; and
\code{points} for the \eqn{(x,y)} data points.

Note that additional parameters (\code{\dots}) are passed only to \code{matplot},
and cannot be used to change the color of plotted points.  To do so,
follow \code{lmrd} by an additional call to the \code{points} function,
e.g. follow \code{lmrd(t3, t4)} by \code{points(t3, t4, col="red")}.
}

%\value{ }

\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}

\examples{
data(airquality)
lmrd(samlmu(airquality$Ozone))
# Tweaking a few graphics parameters makes the graph look better
# (in the author's opinion)
lmrd(samlmu(airquality$Ozone), xaxs="i", yaxs="i", las=1)

# An example that illustrates the sampling variability of L-moments
#
# Generate 50 random samples of size 30 from the Gumbel distribution
# - stored in the rows of matrix mm
mm <- matrix(quagum(runif(1500)), nrow=50)
#
# Compute the first four sample L-moments of each sample
# - stored in the rows of matrix aa
aa <- apply(mm, 1, samlmu)
#
# Plot the L-skewness and L-kurtosis values on an L-moment ratio
# diagram that also shows (only) the population L-moment ratios
# of the Gumbel distribution
lmrd(aa[3,], aa[4,], dist="", twopar="G")
}
\keyword{hplot}
