\name{svcf.sficf}
\alias{svcf.sficf}
\title{
Shapley Value for Causal Forests (SVCF) and Shapley Feature Importance for Causal Forests (SFICF)
}
\description{
Using the estimated CATEs from \code{catecf}, computes (1) SVCF, a local contribution measure that characterizes the directional contribution of each feature to an individual treatment effect; and (2) SFICF, a global feature importance measure that extends SVCF to quantify feature importance at the cohort level. 
}
\usage{
svcf.sficf(tau, X, num.trees = 3000, num.trees.tune = 3000, num.rep = 1,
           tune.oob = TRUE, mtry = NULL, min.node.size = NULL,
           seed = NULL, ...)
}
\arguments{
  \item{tau}{
A numeric vector of estimated CATEs 
of length \code{n}, typically obtained from \code{\link{catecf}}.
}
\item{X}{
A numeric matrix or data frame of covariates with \code{n} rows and \code{p} columns.
Column names are treated as feature names and should be valid R identifiers
(e.g., compatible with \code{make.names}) to ensure seamless integration
with downstream computations.
}
  \item{num.trees}{
Number of trees used to fit the final random forest model for computing
Shapley values. (default: \code{num.trees = 3000})
}
  \item{num.trees.tune}{
Number of trees used during out-of-bag tuning of hyperparameters.
Only used when \code{tune.oob = TRUE}. (default: \code{num.trees.tune = 3000})
}
  \item{num.rep}{
A positive integer specifying the number of repeated random forest fits.
If \code{num.rep = 1}, results are computed from a single fit.
If \code{num.rep > 1}, SVCF and SFICF are averaged across repetitions.
(default: \code{num.rep = 1})
}
\item{tune.oob}{
Logical; if \code{TRUE}, performs out-of-bag (OOB) tuning over
\code{mtry} and \code{min.node.size} using a predefined or user-supplied grid.
If \code{FALSE}, the values supplied via \code{mtry} and
\code{min.node.size} are used directly.
If both are \code{NULL}, the default settings of
\code{ranger::ranger} for continuous outcomes are used.
(default: \code{tune.oob = TRUE})
}
\item{mtry}{
Number of variables randomly selected as candidates at each split.
If \code{tune.oob = TRUE} and \code{mtry = NULL}, a default tuning grid is
constructed based on the number of features \code{p}, typically including
values such as \code{sqrt(p)}, \code{p/10}, \code{p/5}, and \code{p/3}.
If \code{tune.oob = TRUE} and \code{mtry} is supplied as a numeric vector,
the provided values are used directly as tuning candidates.
If \code{tune.oob = FALSE} and \code{mtry = NULL}, the default value used by
\code{ranger::ranger} for continuous outcomes is applied.
}
\item{min.node.size}{
Minimum number of observations in a terminal node.
If \code{tune.oob = TRUE} and \code{min.node.size = NULL}, a default tuning
grid (e.g., \code{5} and \code{10}) is used for OOB tuning.
If \code{tune.oob = TRUE} and \code{min.node.size} is supplied as a numeric
vector, the provided values are used directly as tuning candidates.
If \code{tune.oob = FALSE} and \code{min.node.size = NULL}, the default value
used by \code{ranger::ranger} for continuous outcomes is applied.
}
  \item{seed}{
Optional integer seed for reproducibility. When supported by
\code{ranger::ranger}, the seed is passed directly; otherwise
\code{set.seed()} is used before fitting.
}
  \item{\dots}{
Additional arguments passed to \code{ranger::ranger}.
}
}
\details{
This function first fits a random forest model to the estimated CATEs 
and then computes Shapley values for each feature using
\code{treeshap}. The resulting SVCF provides local, signed feature
contributions for each observation, while SFICF aggregates these contributions
to quantify global feature importance.

When \code{num.rep > 1}, the entire procedure is repeated multiple times using
different random seeds (when provided), and both SVCF and SFICF are averaged
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
A list with the following components:
\item{svcf}{An \code{n x p} numeric matrix of Shapley values representing
local feature contributions for each observation.}
\item{sficf}{A numeric vector of length \code{p} containing global
feature importance scores.}
\item{tau}{The input vector of estimated CATEs.}
\item{res}{A data frame summarizing out-of-bag tuning results when
\code{tune.oob = TRUE}; \code{NULL} otherwise.}
}
\references{
Shapley, L. S. (1953).
\emph{A value for n-person games}.
In \emph{Contributions to the Theory of Games}, Vol.~2, 307--317.
Princeton University Press, Princeton, NJ.

Wager, S. and Athey, S. (2018).
\emph{Estimation and inference of heterogeneous treatment effects using random forests}.
Journal of the American Statistical Association, \bold{113}(523), 1228--1242.

Athey, S., Tibshirani, J., and Wager, S. (2019).
\emph{Generalized random forests}.
The Annals of Statistics, \bold{47}(2), 1148--1178.

Breiman, L. (2001).
\emph{Random forests}.
Machine Learning, \bold{45}, 5--32.

Lundberg, S. M., Erion, G., Chen, H., and others (2020).
\emph{From local explanations to global understanding with explainable AI for trees}.
Nature Machine Intelligence, \bold{2}(1), 56--67.

Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\author{
Hyunwook Koh
}
\note{
This function requires the \pkg{ranger} and \pkg{treeshap} packages.
}
\seealso{
\code{\link{catecf}}, \code{\link[ranger]{ranger}}, \code{\link[treeshap]{treeshap}}
}
\examples{
library(grf)
library(ranger)
library(treeshap)

data(ecigarette)

tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
              num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X, num.rep = 2, seed = 123)

names(svcf.sficf.out)
}
