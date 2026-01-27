# shapcf

Title: Principled Feature Importance and Explanations for Causal Forests via Shapley Values

Version: 1.0

Date: 2026-01-27

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package provides principled, explainable feature attribution methods for causal forests, grounded in the fair allocation principle of the Shapley value. The package implements (1) SVCF (Shapley Values for Causal Forests), a local explanation method that quantifies the signed contribution of each feature to an individual treatment effect, and (2) SFICF (Shapley Feature Importance for Causal Forests), a global feature importance measure obtained by aggregating SVCF across a cohort. In addition, the package provides SHAP-style visualization tools for both local and global explanations, including summary (beeswarm), force (waterfall), and importance (bar) plots.

Depends: R(>= 4.4.1), grf, ranger, treeshap, ggplot2, ggbeeswarm

License: GPL-3

## Reference

* Koh H. Principled Feature Importance and Explanations for Causal Forests via Shapley Values. (_In Review_)

## Troubleshooting Tips

If you have any problems using this R package, please report them via Issues (https://github.com/hk1785/shapcf/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Prerequisites

grf, ranger, treeshap, ggplot2, ggbeeswarm
```
install.packages(c("grf", "ranger", "treeshap", "ggplot2", "ggbeeswarm"))
```

## Installation

```
library(devtools)
install_github("hk1785/shapcf", force = TRUE)
```

---------------------------------------------------------------------------------------------------------------------------------------

# Main Functions

## :mag: catecf

### Description
Estimates conditional average treatment effects (CATEs) using causal forests.

### Syntax
```
catecf(Y, X, W, num.trees = 3000, num.rep = 1, seed = NULL,
tune.parameters = "all", compute.oob.predictions = TRUE, ...)
```

### Arguments
* _Y_ - A numeric vector of outcomes of length _n_.
* _X_ - A numeric matrix or data frame of covariates with _n_ rows and _p_ columns. Column names are treated as feature names and should be valid R identifiers (e.g., compatible with `make.names`).
* _W_ - A treatment indicator vector of length _n_ (typically binary, e.g., 0/1).
* _num.trees_ - Number of trees to grow in each causal forest fit (Default: 3000).
* _num.rep_ - Number of repeated fits. If `num.rep > 1`, returns row-wise average CATEs (Default: 1).
* _seed_ - Optional integer seed for reproducibility.
* _tune.parameters_ - Tuning option passed to `grf::causal_forest` (Default: `"all"`).
* _compute.oob.predictions_ - Logical; whether to compute out-of-bag predictions (Default: `TRUE`).
* _..._ - Additional arguments passed to `grf::causal_forest`.

### Values
A numeric vector of length n containing estimated CATEs for each observation.

### Practical Recommendation
As a practical recommendation, users are advised to first confirm that the function runs correctly with relatively light default settings. Once the workflow is verified, more computationally intensive configurations, such as larger values of num.trees = 30000 and num.rep = 30 with a fixed seed, can be used to obtain more stable and reproducible results.

### Example
Import requisite R packages
```
library(grf)
library(shapcf)
```
Example Data: Oral microbiome data on e-cigarette use and gingival inflammation (see Koh (In Review) for details).
```
data(ecigarette)
```
Estimate CATEs
```
tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
              num.rep = 2, seed = 123)
head(tau)
```
More Details
```
?catecf
```

## :mag: svcf.sficf

### Description
Computes Shapley-based explanations for causal forests using estimated CATEs. Specifically, this function returns (1) **SVCF** (Shapley Values for Causal Forests), which quantify local, signed feature contributions to individual treatment effects, and (2) **SFICF** (Shapley Feature Importance for Causal Forests), a global feature importance measure obtained by aggregating SVCF across individuals.

### Syntax
```
svcf.sficf(tau, X,
num.trees = 3000, num.trees.tune = 3000, num.rep = 1,
tune.oob = TRUE, mtry = NULL, min.node.size = NULL,
seed = NULL, ...)
```

### Arguments
* _tau_ - A numeric vector of estimated CATEs of length _n_, typically obtained from `catecf`.
* _X_ - A numeric matrix or data frame of covariates with _n_ rows and _p_ columns.  
  Column names are treated as feature names and should be valid R identifiers
  (e.g., compatible with `make.names`).
* _num.trees_ - Number of trees used to fit the final random forest model for computing
  Shapley values (Default: 3000).
* _num.trees.tune_ - Number of trees used during out-of-bag tuning of hyperparameters.
  Only used when `tune.oob = TRUE` (Default: 3000).
* _num.rep_ - Number of repeated random forest fits. If `num.rep > 1`, SVCF and SFICF
  are averaged across repetitions (Default: 1).
* _tune.oob_ - Logical; whether to perform out-of-bag tuning over `mtry` and
  `min.node.size`. If `FALSE` and both tuning parameters are `NULL`, the default
  settings of `ranger::ranger` for continuous outcomes are used (Default: TRUE).
* _mtry_ - Number of variables randomly selected as candidates at each split.
  If `tune.oob = TRUE`, this can be a numeric vector specifying tuning candidates.
* _min.node.size_ - Minimum number of observations in a terminal node.
  If `tune.oob = TRUE`, this can be a numeric vector specifying tuning candidates.
* _seed_ - Optional integer seed for reproducibility.
* _..._ - Additional arguments passed to `ranger::ranger`.

### Values
A list with the following components:
* _svcf_ - An _n × p_ numeric matrix of SVCF values representing local feature
  contributions for each observation.
* _sficf_ - A named numeric vector of length _p_ containing normalized global
  feature importance scores.
* _tau_ - The input vector of estimated CATEs.
* _res_ - A data frame summarizing out-of-bag tuning results when `tune.oob = TRUE`
  (otherwise `NULL`).

### Practical Recommendation
As a practical recommendation, users are advised to first confirm that the function runs correctly with relatively light default settings. Once the workflow is verified, more computationally intensive configurations, such as larger values of num.trees = 30000 and num.rep = 30 with a fixed seed, can be used to obtain more stable and reproducible results.

### Example
Import requisite R packages
```
library(grf)
library(ranger)
library(treeshap)
library(shapcf)
```
Example Data: Oral microbiome data on e-cigarette use and gingival inflammation (see Koh (In Review) for details).
```
data(ecigarette)
```
Estimate CATEs and compute SVCF / SFICF
```
tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X,
num.rep = 2, seed = 123)

names(svcf.sficf.out)
```
More Details
```
?svcf.sficf
```

---------------------------------------------------------------------------------------------------------------------------------------

# Visualization Tools

## :mag: beeswarm.svcf.sficf

### Description
Creates a SHAP-style summary (beeswarm) plot for SVCF values, with points colored by the corresponding feature values in `X`.   Features are ordered by `sficf`, and the top `k` features are displayed.

### Syntax
```
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
```

### Arguments
* _svcf_ - An _n × p_ numeric matrix (or data frame) of SVCF values, typically from `svcf.sficf`. Columns should correspond to features.
* _sficf_ - A named numeric vector of length _p_ containing SFICF scores, typically from `svcf.sficf`. Features are ranked by this vector.
* _X_ - A numeric matrix or data frame of covariates with _n_ rows and _p_ columns. Column names must overlap with `colnames(svcf)` and are used for coloring.
* _k_ - Number of top features (ranked by `sficf`) to display (Default: 10).
* _title_ - Optional plot title (Default: `NULL`).
* _xlab_ - Label for the x-axis (Default: `"SVCF"`).
* _point.size_ - Point size for the beeswarm layer (Default: 1.6).
* _point.alpha_ - Point transparency in `[0,1]` (Default: 0.55).
* _low.col_ - Low-end color for feature value gradient (Default: `"#0052A5"`).
* _high.col_ - High-end color for feature value gradient (Default: `"#DC2626"`).
* _qlims_ - Quantile limits in `[0,1]` for truncating the color scale (Default: `c(0.05, 0.95)`).
* _color.alpha_ - Alpha applied to the gradient colors via `adjustcolor()` (Default: 0.60).
* _legend.position_ - Legend position passed to `ggplot2::theme()` (Default: `"bottom"`).
* _legend.direction_ - Legend direction passed to `ggplot2::guide_colorbar()` (Default: `"horizontal"`).
* _legend.barwidth_ - Width of colorbar as a `grid::unit` object (Default: `grid::unit(170, "pt")`).
* _legend.barheight_ - Height of colorbar as a `grid::unit` object (Default: `grid::unit(10, "pt")`).
* _legend.text.size_ - Legend text size (Default: 9).
* _legend.title.size_ - Legend title size (Default: 9).
* _axis.title.x.size_ - X-axis title font size (Default: 11).
* _axis.text.x.size_ - X-axis tick label size (Default: 10).
* _axis.text.y.size_ - Y-axis tick label size (Default: 12).
* _plot.title.size_ - Plot title size (Default: 16).
* _left.margin_ - Left margin to accommodate long feature names (Default: 8).
* _margin.t_ - Top margin (Default: 6).
* _margin.r_ - Right margin (Default: 6).
* _margin.b_ - Bottom margin (Default: 6).

### Values
A `ggplot` object representing the summary (beeswarm) plot.

### Example
Import requisite R packages
```
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)
library(shapcf)
```
Example Data: Oral microbiome data on e-cigarette use and gingival inflammation (see Koh (In Review) for details).
```
data(ecigarette)
```
Estimate CATEs and compute SVCF / SFICF
```
tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X,
num.rep = 2, seed = 123)

names(svcf.sficf.out)
```
Draw beeswarm plot
```
p.beeswarm <- beeswarm.svcf.sficf(svcf = svcf.sficf.out$svcf,
sficf = svcf.sficf.out$sficf,
X = ecigarette$X, k = 20)
p.beeswarm
```

## :mag: waterfall.svcf.sficf

### Description
Creates a SHAP-style force (waterfall) plot that provides a local explanation of feature contributions for a given individual `i`. The plot visualizes signed SVCF values as directional shifts that move the treatment effect from a global baseline (`base.value`) to the individual-specific treatment effect.

### Syntax
```
waterfall.svcf.sficf(svcf, X = NULL, i = 1, k = 10, base.value = 0,
title = NULL, xlab = "Treatment Effect (Cumulative)",
pos.col = "#E07A00", neg.col = "#10978F",
rect.alpha = 0.85, seg.alpha = 0.8, seg.size = 0.9,
rect.h = 0.35, base.linetype = "dashed", base.alpha = 0.5,
digits = 3, value.digits = 5, value.thresh = 1e-05,
plot.title.size = 15, axis.title.x.size = 10,
axis.text.x.size = 10, axis.text.y.size = 12,
left.margin = 8, margin.t = 6, margin.r = 6, margin.b = 6)
```

### Arguments
* _svcf_ - An _n × p_ numeric matrix (or data frame) of SVCF values, typically from `svcf.sficf`. Columns should correspond to features.
* _X_ - Optional feature matrix/data frame with _n_ rows and _p_ columns. If provided, feature values for individual `i` are shown in parentheses in the y-axis labels (excluding `"Other"`).
* _i_ - Index of the observation to explain (Default: 1).
* _k_ - Number of top features (by absolute SVCF) to display; remaining features are aggregated into `"Other"` (Default: 10).
* _base.value_ - Baseline value from which contributions are accumulated. In applications, this is often chosen as the global average treatment effect, e.g., `mean(tau)` (Default: 0).
* _title_ - Optional plot title (Default: `NULL`).
* _xlab_ - Label for the x-axis (Default: `"Treatment Effect (Cumulative)"`).
* _pos.col_ - Fill color for positive contributions (Default: `"#E07A00"`).
* _neg.col_ - Fill color for negative contributions (Default: `"#10978F"`).
* _rect.alpha_ - Alpha for contribution rectangles in `[0,1]` (Default: 0.85).
* _seg.alpha_ - Alpha for connecting segments in `[0,1]` (Default: 0.8).
* _seg.size_ - Line width for connecting segments (Default: 0.9).
* _rect.h_ - Half-height of the contribution rectangles (Default: 0.35).
* _base.linetype_ - Line type for baseline reference line (Default: `"dashed"`).
* _base.alpha_ - Alpha for baseline reference line in `[0,1]` (Default: 0.5).
* _digits_ - Digits for x-axis tick labels (Default: 3).
* _value.digits_ - Digits for feature values in labels when `X` is provided (Default: 5).
* _value.thresh_ - Threshold below which feature values are shown as `<value.thresh` (Default: 1e-05).
* _plot.title.size_ - Plot title size (Default: 15).
* _axis.title.x.size_ - X-axis title size (Default: 10).
* _axis.text.x.size_ - X-axis tick label size (Default: 10).
* _axis.text.y.size_ - Y-axis tick label size (Default: 12).
* _left.margin_ - Left margin to accommodate long labels (Default: 8).
* _margin.t_ - Top margin (Default: 6).
* _margin.r_ - Right margin (Default: 6).
* _margin.b_ - Bottom margin (Default: 6).

### Values
A `ggplot` object representing the force (waterfall) plot.

### Example
Import requisite R packages
```
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)
library(shapcf)
```
Example Data: Oral microbiome data on e-cigarette use and gingival inflammation (see Koh (In Review) for details).
```
data(ecigarette)
```
Estimate CATEs and compute SVCF / SFICF
```
tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X,
num.rep = 2, seed = 123)

names(svcf.sficf.out)
```
Draw waterfall plot (choose baseline as mean CATE)
```
p.waterfall <- waterfall.svcf.sficf(svcf = svcf.sficf.out$svcf,
X = ecigarette$X,
i = 1, k = 10,
base.value = mean(tau))
p.waterfall
```

## :mag: bar.sficf

### Description
Creates a horizontal bar plot of the top `k` features ranked by SFICF (global feature importance) scores.

### Syntax
```
bar.sficf(sficf, k = 10, title = NULL, xlab = "Feature Importance",
bar.col = "#6B7280", bar.alpha = 0.9, bar.width = 0.7,
base.size = 12, plot.title.size = 12,
axis.title.x.size = 11, axis.text.x.size = 10, axis.text.y.size = 11,
left.margin = 8, margin.t = 6, margin.r = 6, margin.b = 6)
```

### Arguments
* _sficf_ - A named numeric vector of length _p_ containing SFICF scores, typically from `svcf.sficf`.
* _k_ - Number of top features to display (Default: 10).
* _title_ - Optional plot title (Default: `NULL`).
* _xlab_ - Label for the x-axis (Default: `"Feature Importance"`).
* _bar.col_ - Fill color for bars (Default: `"#6B7280"`).
* _bar.alpha_ - Bar transparency in `[0,1]` (Default: 0.9).
* _bar.width_ - Bar width passed to `ggplot2::geom_col` (Default: 0.7).
* _base.size_ - Base font size for `ggplot2::theme_bw` (Default: 12).
* _plot.title.size_ - Plot title size (Default: 12).
* _axis.title.x.size_ - X-axis title size (Default: 11).
* _axis.text.x.size_ - X-axis tick label size (Default: 10).
* _axis.text.y.size_ - Y-axis tick label size (Default: 11).
* _left.margin_ - Left margin to accommodate long labels (Default: 8).
* _margin.t_ - Top margin (Default: 6).
* _margin.r_ - Right margin (Default: 6).
* _margin.b_ - Bottom margin (Default: 6).

### Values
A `ggplot` object representing the importance (bar) plot.

### Example
Import requisite R packages
```
library(grf)
library(ranger)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)
library(shapcf)
```
Example Data: Oral microbiome data on e-cigarette use and gingival inflammation (see Koh (In Review) for details).
```
data(ecigarette)
```
Estimate CATEs and compute SVCF / SFICF
```
tau <- catecf(Y = ecigarette$Y, X = ecigarette$X, W = ecigarette$W,
num.rep = 2, seed = 123)

svcf.sficf.out <- svcf.sficf(tau = tau, X = ecigarette$X,
num.rep = 2, seed = 123)

names(svcf.sficf.out)
```
Draw bar plot
```
p.bar <- bar.sficf(sficf = svcf.sficf.out$sficf, k = 20)
p.bar
```
---------------------------------------------------------------------------------------------------------------------------------------

# Example Datasets

## :mag: ecigarette

### Description
Subgingival microbiome data used to study how microbial genera moderate the effect of e-cigarette use on gingival inflammation among participants from Baltimore, MD.

### Usage
```
data(ecigarette)
```

### Format
A list with three components:
* _Y_ - Binary outcome indicating gingival inflammation (0 = healthy, 1 = diseased).
* _W_ - Binary exposure indicator (0 = non-use, 1 = e-cigarette use).
* _X_ - A data frame of relative abundances for `p = 71` subgingival microbial genera and `n = 197` participants.

### References
Park, B., Koh, H., Patatanian, M., and others (2023). _The mediating roles of the oral microbiome in saliva and subgingival sites between e-cigarette smoking and gingival inflammation_. BMC Microbiology, 23, 35. doi:10.1186/s12866-023-02779-z.

Koh H. Principled Feature Importance and Explanations for Causal Forests via Shapley Values. (_In Review_)

## :mag: antibiotic

### Description
Gut microbiome data from a mouse study investigating how antibiotic treatment (tylosin; 0 = control, 1 = treatment) moderates the onset of type 1 diabetes (0 = healthy, 1 = diseased).

### Usage
```
data(antibiotic)
```

### Format
A list with three components:
* _Y_ - Binary outcome indicating type 1 diabetes onset (0 = healthy, 1 = diseased).
* _W_ - Binary treatment indicator for antibiotic exposure (0 = control, 1 = treatment).
* _X_ - A data frame of relative abundances for `p = 29` gut microbial genera and `n = 521` observations.

### References
Zhang, X. S., Li, J., Krautkramer, K. A., et al. (2018). _Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity_. eLife, 7, e37816.

Koh H. Principled Feature Importance and Explanations for Causal Forests via Shapley Values. (_In Review_)

## :mag: immuno

### Description
Gut microbiome data from metastatic melanoma patients, used to study how microbial genera moderate the effect of cancer immunotherapy on clinical response.

### Usage
```
data(immuno)
```

### Format
A list with three components:
* _Y_ - Binary outcome indicating clinical response (0 = non-response, 1 = response).
* _W_ - Binary treatment indicator (0 = anti--PD-1 monotherapy, 1 = combined anti--PD-1/anti--CTLA-4 therapy).
* _X_ - A data frame of relative abundances for `p = 119` gut microbial genera and `n = 255` patients.

### References
Limeta, A., Ji, B., Levin, M., Gatto, F., and Nielsen, J. (2020). _Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma_. JCI Insight, 5(23), e140940.

Koh H. Principled Feature Importance and Explanations for Causal Forests via Shapley Values. (_In Review_)













