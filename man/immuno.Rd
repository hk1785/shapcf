\name{immuno}
\alias{immuno}
\docType{data}
\title{
Cancer immunotherapy response microbiome data
}
\description{
Gut microbiome data from metastatic melanoma patients, used to study how microbial
genera moderate the effect of cancer immunotherapy on clinical response.
}
\usage{data("immuno")}
\format{
A list with three components:
\describe{
  \item{Y}{Binary outcome indicating clinical response (0 = non-response, 1 = response).}
  \item{W}{Binary treatment indicator (0 = anti--PD-1 monotherapy, 1 = combined anti--PD-1/anti--CTLA-4 therapy).}
  \item{X}{A data frame of relative abundances for \code{p = 119} gut microbial genera.}
}
The dataset contains \code{n = 255} patients.
}
\details{
The data integrate multiple studies, with study-specific variation adjusted via
ConQuR prior to analysis. Further preprocessing and analysis details are provided
in Koh (In Review).
}
\source{
Originally reported in Limeta et al. (2020).
}
\references{
Limeta, A., Ji, B., Levin, M., Gatto, F., and Nielsen, J. (2020).
\emph{Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma}.
JCI Insight, \bold{5}(23), e140940.

Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\examples{
data(immuno)
str(immuno)
}
\keyword{datasets}