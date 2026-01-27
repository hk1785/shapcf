\name{ecigarette}
\alias{ecigarette}
\docType{data}
\title{
E-cigarette use and gingival inflammation microbiome data
}
\description{
Subgingival microbiome data used to study how microbial genera moderate the effect
of e-cigarette use on gingival inflammation among participants from Baltimore, MD.
}
\usage{data("ecigarette")}
\format{
A list with three components:
\describe{
  \item{Y}{Binary outcome indicating gingival inflammation (0 = healthy, 1 = diseased).}
  \item{W}{Binary exposure indicator (0 = non-use, 1 = e-cigarette use).}
  \item{X}{A data frame of relative abundances for \code{p = 71} subgingival microbial genera.}
}
The dataset contains \code{n = 197} participants.
}
\details{
The data include two sampling time points, with time-specific variation adjusted
via ConQuR prior to analysis. Further preprocessing and analysis details are
provided in Koh (In Review).
}
\source{
Originally reported in Park et al. (2023).
}
\references{
Park, B., Koh, H., Patatanian, M., and others (2023).
\emph{The mediating roles of the oral microbiome in saliva and subgingival sites between e-cigarette smoking and gingival inflammation}.
BMC Microbiology, \bold{23}, 35. doi:10.1186/s12866-023-02779-z.

Koh, H. (In Review).
\emph{Principled Feature Importance and Explanations for Causal Forests via Shapley Values}.
}
\examples{
data(ecigarette)
str(ecigarette)
}
\keyword{datasets}