\name{catecf}
\alias{catecf}
\title{
CATE estimation via causal forests
}
\description{
Estimates conditional average treatment effects (CATEs) using \code{grf::causal_forest}. Optionally repeats model fitting multiple times and returns 
the averaged CATEs to reduce Monte Carlo variability.
}
\usage{
catecf(Y, X, W, num.trees = 3000, num.rep = 1, seed = NULL,
       tune.parameters = "all", compute.oob.predictions = TRUE, ...)
}
\arguments{
  \item{Y}{
A numeric vector of outcomes of length \code{n}.
}
\item{X}{
A numeric matrix or data frame of covariates with \code{n} rows and \code{p} columns.
Column names are treated as feature names and should be valid R identifiers
(e.g., compatible with \code{make.names}) to ensure seamless integration
with downstream computations.
}
  \item{W}{
A treatment indicator vector of length \code{n} (typically binary, e.g., 0/1).
}
  \item{num.trees}{
Number of trees to grow in each causal forest fit. (default: \code{num.trees = 3000})
}
  \item{num.rep}{
A positive integer. If \code{num.rep = 1}, fits a single causal forest.
If \code{num.rep > 1}, fits \code{num.rep} forests and returns the row-wise average
of predicted CATEs. (default: \code{num.rep = 1})
}
  \item{seed}{
Optional integer seed for reproducibility. If the installed \code{grf::causal_forest}
supports a \code{seed} argument, it is passed directly; otherwise \code{set.seed()}
is used before fitting.
}
  \item{tune.parameters}{
Tuning option passed to \code{grf::causal_forest} (e.g., \code{"all"}, \code{"none"}).
See \code{\link[grf]{causal_forest}} for details.
}
  \item{compute.oob.predictions}{
Logical; passed to \code{grf::causal_forest}. If \code{TRUE}, out-of-bag predictions
are computed. (default: \code{compute.oob.predictions = TRUE})
}
  \item{\dots}{
Additional arguments passed to \code{grf::causal_forest}.
}
}
\details{
When \code{num.rep > 1}, the entire procedure is repeated multiple times using
different random seeds (when provided), and CATEs are averaged
across repetitions. This can help stabilize the results by reducing Monte Carlo
variability.

As a practical recommendation, users are advised to first confirm that the
function runs correctly with relatively light default settings. Once the
workflow is verified, more computationally intensive configurations, such as
larger values of \code{num.trees = 30000} and \code{num.rep = 50}
with a fixed \code{seed}, can be used to obtain more stable and reproducible
results.
}
\value{
A numeric vector of length \code{n} containing estimated CATEs for each observation.
}
\references{
Wager, S. and Athey, S. (2018).
\emph{Estimation and inference of heterogeneous treatment effects using random forests}.
Journal of the American Statistical Association, \bold{113}(523), 1228--1242.

Athey, S., Tibshirani, J., and Wager, S. (2019).
\emph{Generalized random forests}.
The Annals of Statistics, \bold{47}(2), 1148--1178.

Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\author{
Hyunwook Koh
}
\note{
This function requires the \pkg{grf} package.
}
\seealso{
\code{\link[grf]{causal_forest}}, \code{\link[grf]{predict.causal_forest}}
}
\examples{
library(grf)

data(ecigarette)

tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W, num.rep = 2, seed = 123)
tau
}