\name{antibiotic}
\alias{antibiotic}
\docType{data}
\title{
Antibiotic treatment and type 1 diabetes microbiome data
}
\description{
Gut microbiome data from a mouse study investigating how antibiotic treatment
(tylosin; 0 = control, 1 = treatment) moderates the onset of type 1 diabetes
(0 = healthy, 1 = diseased).
}
\usage{data("antibiotic")}
\format{
A list with three components:
\describe{
  \item{Y}{Binary outcome indicating type 1 diabetes onset (0 = healthy, 1 = diseased).}
  \item{W}{Binary treatment indicator for antibiotic exposure (0 = control, 1 = treatment).}
  \item{X}{A data frame of relative abundances for \code{p = 29} gut microbial genera.}
}
The dataset contains \code{n = 521} observations.
}
\details{
The data were collected from non-obese diabetic mice across multiple sampling
time points and both sexes. Time- and sex-specific effects were adjusted using
ConQuR prior to analysis. Further details on preprocessing and analysis are
provided in Koh (In Review).
}
\source{
Originally reported in Zhang et al. (2018).
}
\references{
Zhang, X. S., Li, J., Krautkramer, K. A., et al. (2018).
\emph{Antibiotic-induced acceleration of type 1 diabetes alters maturation of
innate intestinal immunity}.
eLife, \bold{7}, e37816.

Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via
Shapley Values}.
}
\examples{
data(antibiotic)
str(antibiotic)
}
\keyword{datasets}