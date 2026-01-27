\name{shapcf-package}
\alias{shapcf-package}
\alias{shapcf}
\docType{package}
\title{
\packageTitle{shapcf}
}
\description{
\packageDescription{shapcf}
}
\details{
The \pkg{shapcf} package provides Shapley-value-based explanation tools for causal
forests. The main workflow is:
(1) estimate CATEs with \code{\link{catecf}},
(2) compute SVCF/SFICF with \code{\link{svcf.sficf}},
and (3) visualize results with \code{\link{beeswarm.svcf.sficf}},
\code{\link{waterfall.svcf.sficf}}, and \code{\link{bar.sficf}}.
}
\author{
\packageAuthor{shapcf}

Maintainer: \packageMaintainer{shapcf}
}
\references{
Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\seealso{
\code{\link{catecf}}, \code{\link{svcf.sficf}},
\code{\link{beeswarm.svcf.sficf}}, \code{\link{waterfall.svcf.sficf}}, \code{\link{bar.sficf}}
}
\keyword{package}