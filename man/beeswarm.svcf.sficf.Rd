\name{beeswarm.svcf.sficf}
\alias{beeswarm.svcf.sficf}
\title{
Summary (beeswarm) plot for SVCF and SFICF
}
\description{
Creates a SHAP-style summary (beeswarm) plot for \code{svcf} values, with points
colored by the corresponding feature values in \code{X}. Features are ordered by
\code{sficf} and the top \code{k} features are displayed.
}
\usage{
beeswarm.svcf.sficf(svcf, sficf, X, k = 10,
                    title = NULL, xlab = "SVCF",
                    point.size = 1.6, point.alpha = 0.55,
                    low.col = "#0052A5", high.col = "#DC2626",
                    qlims = c(0.05, 0.95), color.alpha = 0.60,
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.barwidth = grid::unit(170, "pt"),
                    legend.barheight = grid::unit(10, "pt"),
                    legend.text.size = 9,
                    legend.title.size = 9,
                    axis.title.x.size = 11,
                    axis.text.x.size = 10,
                    axis.text.y.size = 12,
                    plot.title.size = 16,
                    left.margin = 8,
                    margin.t = 6, margin.r = 6, margin.b = 6)
}
\arguments{
  \item{svcf}{
An \code{n x p} numeric matrix (or data frame) of SVCF, typically
from \code{\link{svcf.sficf}}. Columns should correspond to features.
}
\item{sficf}{
A named numeric vector of length \code{p} containing SFICF scores, typically
obtained from \code{\link{svcf.sficf}}. Features are ranked by this vector to
determine the top \code{k} features.
}
\item{X}{
A numeric matrix or data frame of covariates with \code{n} rows and \code{p} columns.
Column names are treated as feature names and should be valid R identifiers
(e.g., compatible with \code{make.names}) to ensure seamless integration
with downstream computations. Column names must overlap with
\code{colnames(svcf)} and are used to color the points in the plot.
}
  \item{k}{
A positive integer specifying the number of top features (ranked by \code{sficf})
to display. (default: \code{k = 10})
}
  \item{title}{
Optional plot title. If \code{NULL}, no title is displayed. (default: \code{title = NULL})
}
  \item{xlab}{
Label for the x-axis. (default: \code{xlab = "SVCF"})
}
  \item{point.size}{
Point size for the beeswarm layer. (default: \code{point.size = 1.6})
}
  \item{point.alpha}{
Transparency of points in the beeswarm layer. Must be between 0 and 1.
(default: \code{point.alpha = 0.55})
}
  \item{low.col}{
Low-end color used for the feature value gradient. (default: \code{low.col = "#0052A5"})
}
  \item{high.col}{
High-end color used for the feature value gradient. (default: \code{high.col = "#DC2626"})
}
  \item{qlims}{
A length-2 numeric vector in \code{[0,1]} specifying quantile limits for the color
scale, used to reduce the influence of extreme feature values. Values outside
these limits are squished. (default: \code{qlims = c(0.05, 0.95)})
}
  \item{color.alpha}{
Alpha used when applying \code{adjustcolor()} to \code{low.col} and \code{high.col}.
Must be between 0 and 1. (default: \code{color.alpha = 0.60})
}
  \item{legend.position}{
Legend position passed to \code{ggplot2::theme()}, e.g., \code{"bottom"}, \code{"right"},
or \code{"none"}. (default: \code{legend.position = "bottom"})
}
  \item{legend.direction}{
Legend direction passed to \code{ggplot2::guide_colorbar()}, typically
\code{"horizontal"} or \code{"vertical"}. (default: \code{legend.direction = "horizontal"})
}
  \item{legend.barwidth}{
Width of the colorbar, specified as a \code{grid::unit} object.
(default: \code{grid::unit(170, "pt")})
}
  \item{legend.barheight}{
Height of the colorbar, specified as a \code{grid::unit} object.
(default: \code{grid::unit(10, "pt")})
}
  \item{legend.text.size}{
Font size for legend text. (default: \code{legend.text.size = 9})
}
  \item{legend.title.size}{
Font size for legend title. (default: \code{legend.title.size = 9})
}
  \item{axis.title.x.size}{
Font size for the x-axis title. (default: \code{axis.title.x.size = 11})
}
  \item{axis.text.x.size}{
Font size for x-axis tick labels. (default: \code{axis.text.x.size = 10})
}
  \item{axis.text.y.size}{
Font size for y-axis tick labels (feature names). (default: \code{axis.text.y.size = 12})
}
  \item{plot.title.size}{
Font size for the plot title. (default: \code{plot.title.size = 16})
}
  \item{left.margin}{
Left plot margin (in \code{ggplot2::margin} units) to provide space for long feature names.
(default: \code{left.margin = 8})
}
  \item{margin.t}{
Top plot margin. (default: \code{margin.t = 6})
}
  \item{margin.r}{
Right plot margin. (default: \code{margin.r = 6})
}
  \item{margin.b}{
Bottom plot margin. (default: \code{margin.b = 6})
}
}
\details{
The function aligns features by taking the intersection of \code{colnames(svcf)}
and \code{colnames(X)}, then selects the top \code{k} features according to
\code{sficf} (using the overlap between \code{names(sficf)} and \code{colnames(svcf)}).
The beeswarm plot is drawn using \code{ggbeeswarm::geom_quasirandom},
and point colors represent feature values in \code{X} with a truncated (quantile-based)
color scale controlled by \code{qlims}.
}
\value{
A \code{ggplot} object representing the summary (beeswarm) plot.
}
\references{
Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\author{
Hyunwook Koh
}
\note{
This function requires the \pkg{ggplot2} and \pkg{ggbeeswarm} packages.
}
\seealso{
\code{\link{svcf.sficf}}, \code{\link[ggbeeswarm]{geom_quasirandom}}, \code{\link[ggplot2]{ggplot}}
}
\examples{
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)

data(ecigarette)

tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
              num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X, num.rep = 2, seed = 123)

beeswarm.out <- beeswarm.svcf.sficf(svcf = svcf.sficf.out$svcf,
                                    sficf = svcf.sficf.out$sficf,
                                    X = ecigarette$X, k = 20)
beeswarm.out
}
