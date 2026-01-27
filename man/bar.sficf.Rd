\name{bar.sficf}
\alias{bar.sficf}
\title{
Importance (bar) plot for SFICF
}
\description{
Creates a horizontal bar plot of the top \code{k} features ranked by SFICF
(global feature importance) scores.
}
\usage{
bar.sficf(sficf, k = 10, title = NULL, xlab = "Feature Importance",
          bar.col = "#6B7280", bar.alpha = 0.9, bar.width = 0.7,
          base.size = 12, plot.title.size = 12,
          axis.title.x.size = 11, axis.text.x.size = 10, axis.text.y.size = 11,
          left.margin = 8, margin.t = 6, margin.r = 6, margin.b = 6)
}
\arguments{
  \item{sficf}{
A named numeric vector of length \code{p} containing SFICF scores, typically
obtained from \code{\link{svcf.sficf}}. Names are treated as feature names and
should be valid R identifiers (e.g., compatible with \code{make.names}).
}
  \item{k}{
A positive integer specifying the number of top features to display.
(default: \code{k = 10})
}
  \item{title}{
Optional plot title. If \code{NULL}, no title is displayed. (default: \code{title = NULL})
}
  \item{xlab}{
Label for the x-axis. (default: \code{xlab = "Feature Importance"})
}
  \item{bar.col}{
Fill color for bars. (default: \code{bar.col = "#6B7280"})
}
  \item{bar.alpha}{
Transparency of bars. Must be between 0 and 1. (default: \code{bar.alpha = 0.9})
}
  \item{bar.width}{
Bar width passed to \code{ggplot2::geom_col}. (default: \code{bar.width = 0.7})
}
  \item{base.size}{
Base font size for \code{ggplot2::theme_bw}. (default: \code{base.size = 12})
}
  \item{plot.title.size}{
Font size for the plot title. (default: \code{plot.title.size = 12})
}
  \item{axis.title.x.size}{
Font size for the x-axis title. (default: \code{axis.title.x.size = 11})
}
  \item{axis.text.x.size}{
Font size for x-axis tick labels. (default: \code{axis.text.x.size = 10})
}
  \item{axis.text.y.size}{
Font size for y-axis tick labels (feature names). (default: \code{axis.text.y.size = 11})
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
The function sorts \code{sficf} in decreasing order, selects the top \code{k}
features, and draws a horizontal bar plot using \code{ggplot2::geom_col}.
}
\value{
A \code{ggplot} object representing the importance (bar) plot.
}
\references{
Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\author{
Hyunwook Koh
}
\note{
This function requires the \pkg{ggplot2} package.
}
\seealso{
\code{\link{svcf.sficf}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_col}}
}
\examples{
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)

data(immuno)

tau <- catecf(Y = immuno$Y, X = immuno$X, W = immuno$W,
              num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = immuno$X, num.rep = 2, seed = 123)

bar.out <- bar.sficf(sficf = svcf.sficf.out$sficf, k = 20)
bar.out
}