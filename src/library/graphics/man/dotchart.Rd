% File src/library/graphics/man/dotchart.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2022 R Core Team
% Distributed under GPL 2 or later

\name{dotchart}
\alias{dotchart}
\title{Cleveland's Dot Plots}
\description{
  Draw a Cleveland dot plot.
}
\usage{
dotchart(x, labels = NULL, groups = NULL, gdata = NULL, offset = 1/8,
         ann = par("ann"), xaxt = par("xaxt"), frame.plot = TRUE, log = "",
         cex = par("cex"), pt.cex = cex,
         pch = 21, gpch = 21, bg = par("bg"),
         color = par("fg"), gcolor = par("fg"), lcolor = "gray",
         xlim = range(x[is.finite(x)]),
         main = NULL, xlab = NULL, ylab = NULL, \dots)
}
\arguments{
  \item{x}{either a vector or matrix of numeric values (\code{NA}s are
    allowed).  If \code{x} is a matrix the overall plot consists of
    juxtaposed dotplots for each row.  Inputs which satisfy
    \code{\link{is.numeric}(x)} but not
    \code{is.vector(x) || is.matrix(x)} are coerced by
    \code{\link{as.numeric}}, with a warning.}
  \item{labels}{a vector of labels for each point.
    For vectors the default is to use \code{names(x)} and for matrices
    the row labels \code{dimnames(x)[[1]]}.}
  \item{groups}{an optional factor indicating how the elements of
    \code{x} are grouped.
    If \code{x} is a matrix, \code{groups} will default to the columns
    of \code{x}.}
  \item{gdata}{data values for the groups.
    This is typically a summary such as the median or mean
    of each group.}
  \item{offset}{offset in inches of \code{ylab} and \code{labels}.}
  \item{ann}{a \code{\link{logical}} value indicating whether the default annotation
    (title and x and y axis labels) should appear on the plot.}
  \item{xaxt}{a string indicating the x-axis style; use \code{"n"} to
    suppress and see also \code{\link{par}("xaxt")}.}
  \item{frame.plot}{a logical indicating whether a box should be drawn
    around the plot.}
  \item{log}{a character string indicating if one or the other axis should
    be logarithmic, see \code{\link{plot.default}}.}
  \item{cex}{the character size to be used.  Setting \code{cex}
    to a value smaller than one can be a useful way of avoiding label
    overlap.  Unlike many other graphics functions, this sets the actual
    size, not a multiple of \code{par("cex")}.}
  \item{pt.cex}{the \code{cex} to be applied to plotting symbols.
    This behaves like \code{cex} in \code{plot()}.}
  \item{pch}{the plotting character or symbol to be used.}
  \item{gpch}{the plotting character or symbol to be used for group
    values.}
  \item{bg}{the background color of plotting characters or symbols to be
    used; use \code{\link{par}(bg= *)} to set the background color of
    the whole plot.}
  \item{color}{the color(s) to be used for points and labels.}
  \item{gcolor}{the single color to be used for group labels and
    values.}
  \item{lcolor}{the color(s) to be used for the horizontal lines.}
  \item{xlim}{horizontal range for the plot, see
    \code{\link{plot.window}}, for example.}
  \item{main}{overall title for the plot, see \code{\link{title}}.}
  \item{xlab, ylab}{axis annotations as in \code{title}.}
  \item{\dots}{\link{graphical parameters} can also be specified as arguments.}
}
\value{
  This function is invoked for its side effect, which is to produce two
  variants of dotplots as described in Cleveland (1985).

  Dot plots are a reasonable substitute for bar plots.
}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.

  Cleveland, W. S. (1985)
  \emph{The Elements of Graphing Data.}
  Monterey, CA: Wadsworth.

  Murrell, P. (2005) \emph{R Graphics}. Chapman & Hall/CRC Press.
}
\examples{
dotchart(VADeaths, main = "Death Rates in Virginia - 1940")

op <- par(xaxs = "i")  # 0 -- 100\%
dotchart(t(VADeaths), xlim = c(0,100), bg = "skyblue",
         main = "Death Rates in Virginia - 1940", xlab = "rate [ \% ]",
         ylab = "Grouping:  Age  x   Urbanity . Gender")
par(op)
}
\keyword{hplot}
