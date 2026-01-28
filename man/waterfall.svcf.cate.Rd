\name{waterfall.svcf.cate}
\alias{waterfall.svcf.cate}
\title{
Force (waterfall) plot for SVCF
}
\description{
Creates a SHAP-style force (waterfall) plot that provides a local explanation of
feature contributions for a given individual \code{i}. The plot visualizes signed
SVCF values as directional shifts that move the treatment effect from a global
baseline to the individual-specific treatment effect.
}
\usage{
waterfall.svcf.cate(svcf, tau, X = NULL, i = 1, k = 10, 
                    base.value = NULL, title = NULL, 
					xlab = "Treatment Effect (Cumulative)", 
                    pos.col = "#E07A00", neg.col = "#10978F", rect.alpha = 0.85, seg.alpha = 0.8, 
                    seg.size = 0.9, rect.h = 0.35, base.linetype = "dashed", 
                    base.alpha = 0.5, digits = 3, value.digits = 5, value.thresh = 1e-05, 
                    plot.title.size = 15, axis.title.x.size = 10, axis.text.x.size = 10, 
                    axis.text.y.size = 12, left.margin = 8, margin.t = 6, margin.r = 6, 
					margin.b = 6)
}
\arguments{
  \item{svcf}{
An \code{n x p} numeric matrix (or data frame) of SVCF, typically
from \code{\link{svcf.sficf}}. Columns should correspond to features.
}
  \item{tau}{
A numeric vector of estimated CATEs 
of length \code{n}, typically obtained from \code{\link{catecf}}.
}
  \item{X}{
Optional feature matrix/data frame with \code{n} rows and \code{p} columns.
If provided, the observed feature values for individual \code{i} are displayed in
parentheses in the y-axis labels (excluding \code{"Other"}), thereby describing the
individual's profile. Column names are treated as feature names and should be valid
R identifiers (e.g., compatible with \code{make.names}) to ensure seamless integration
with downstream computations. Column names must overlap with \code{colnames(svcf)}.
}
  \item{i}{
Index of the observation to explain. Must be between \code{1} and \code{n}.
(default: \code{i = 1})
}
  \item{k}{
A positive integer specifying the number of top features (by absolute SVCF)
to display. Remaining features are aggregated into \code{"Other"}.
(default: \code{k = 10})
}
\item{base.value}{
Baseline value from which feature contributions are accumulated. In practice, this is
typically chosen as the global average treatment effect, for example
\eqn{\frac{1}{n}\sum_{i=1}^n \hat{\tau}(x_i)}.
The default setting is \code{base.value = NULL}, in which case
\code{base.value} is set to \code{mean(tau)}.
}
  \item{title}{
Optional plot title. If \code{NULL}, no title is displayed. (default: \code{title = NULL})
}
  \item{xlab}{
Label for the x-axis. (default: \code{xlab = "Treatment Effect (Cumulative)"})
}
  \item{pos.col}{
Fill color for positive contributions. (default: \code{pos.col = "#E07A00"})
}
  \item{neg.col}{
Fill color for negative contributions. (default: \code{neg.col = "#10978F"})
}
  \item{rect.alpha}{
Alpha level for the contribution rectangles. Must be between 0 and 1.
(default: \code{rect.alpha = 0.85})
}
  \item{seg.alpha}{
Alpha level for the connecting segments. Must be between 0 and 1.
(default: \code{seg.alpha = 0.8})
}
  \item{seg.size}{
Line width for the connecting segments. (default: \code{seg.size = 0.9})
}
  \item{rect.h}{
Half-height of the contribution rectangles. (default: \code{rect.h = 0.35})
}
  \item{base.linetype}{
Line type for the baseline reference line. (default: \code{base.linetype = "dashed"})
}
  \item{base.alpha}{
Alpha level for the baseline reference line. Must be between 0 and 1.
(default: \code{base.alpha = 0.5})
}
  \item{digits}{
Number of digits for x-axis tick labels. (default: \code{digits = 3})
}
  \item{value.digits}{
Number of digits for feature values shown in labels when \code{X} is provided.
(default: \code{value.digits = 5})
}
  \item{value.thresh}{
Threshold below which feature values are displayed as \code{<value.thresh}.
(default: \code{value.thresh = 1e-05})
}
  \item{plot.title.size}{
Font size for the plot title. (default: \code{plot.title.size = 15})
}
  \item{axis.title.x.size}{
Font size for the x-axis title. (default: \code{axis.title.x.size = 10})
}
  \item{axis.text.x.size}{
Font size for x-axis tick labels. (default: \code{axis.text.x.size = 10})
}
  \item{axis.text.y.size}{
Font size for y-axis tick labels (feature names). (default: \code{axis.text.y.size = 12})
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
The function extracts the SVCF vector for observation \code{i} and orders features
by the magnitude of their contributions. The top \code{k} features are displayed
explicitly, while the remaining features are aggregated into \code{"Other"}.
Each bar represents a feature, where bar length indicates contribution magnitude and
bar direction (color) indicates contribution sign. By sequentially aggregating
feature-specific contributions, the plot illustrates how \code{base.value} is
transformed into the individual-specific treatment effect; the cumulative path is
read downward in sequence. If \code{X} is provided, the observed feature values for
individual \code{i} are appended in parentheses to the corresponding feature labels.
}
\value{
A \code{ggplot} object representing the force (waterfall) plot.
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
\code{\link{svcf.sficf}}, \code{\link[ggplot2]{ggplot}}
}
\examples{
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)

data(ecigarette)

tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
              num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X, num.rep = 2, seed = 123)

waterfall.out <- waterfall.svcf.cate(svcf = svcf.sficf.out$svcf, tau = tau,
                                     X = ecigarette$X, i = 1, k = 20)
waterfall.out
}

