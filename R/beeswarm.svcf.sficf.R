beeswarm.svcf.sficf <-
function(svcf, sficf, X, k = 10,
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
                                margin.t = 6, margin.r = 6, margin.b = 6) {
  
  ## ---- checks
  stopifnot(is.matrix(svcf) || is.data.frame(svcf))
  stopifnot(is.matrix(X)    || is.data.frame(X))
  stopifnot(is.numeric(sficf), !is.null(names(sficf)))
  stopifnot(k >= 1, k == as.integer(k))
  stopifnot(length(qlims) == 2, is.numeric(qlims), all(is.finite(qlims)),
            all(qlims >= 0), all(qlims <= 1), qlims[1] < qlims[2])
  stopifnot(point.alpha >= 0, point.alpha <= 1, color.alpha >= 0, color.alpha <= 1)
  
  if (is.null(colnames(svcf)) || is.null(colnames(X))) {
    stop("Both 'svcf' and 'X' must have column names.")
  }
  
  ## ---- align columns
  common <- intersect(colnames(svcf), colnames(X))
  if (length(common) == 0) {
    stop("No common feature names found between 'svcf' and 'X'.")
  }
  svcf <- svcf[, common, drop = FALSE]
  X    <- X[,    common, drop = FALSE]
  
  ## ---- choose top-k vars by sficf
  sficf2 <- sort(sficf, decreasing = TRUE)
  vars <- intersect(names(sficf2), colnames(svcf))
  if (length(vars) == 0) {
    stop("No overlap between names(sficf) and colnames(svcf)/colnames(X).")
  }
  vars <- vars[seq_len(min(k, length(vars)))]
  
  ## ---- build long data
  sv <- as.matrix(svcf[, vars, drop = FALSE])
  xv <- as.matrix(X[,    vars, drop = FALSE])
  
  df <- data.frame(
    id = rep(seq_len(nrow(sv)), times = length(vars)),
    feature = rep(vars, each = nrow(sv)),
    shap = as.numeric(c(sv)),
    value = as.numeric(c(xv)),
    stringsAsFactors = FALSE
  )
  df$feature <- factor(df$feature, levels = rev(vars))
  
  ## ---- color limits
  lims <- stats::quantile(df$value, probs = qlims, na.rm = TRUE, names = FALSE, type = 7)
  if (!all(is.finite(lims))) lims <- c(0, 1)
  
  ggplot2::ggplot(df, ggplot2::aes(x = shap, y = feature, color = value)) +
    ggbeeswarm::geom_quasirandom(groupOnX = TRUE,
                                 size = point.size, alpha = point.alpha) +
    ggplot2::scale_color_gradient(
      low  = grDevices::adjustcolor(low.col,  alpha.f = color.alpha),
      high = grDevices::adjustcolor(high.col, alpha.f = color.alpha),
      limits = lims,
      oob = scales::squish
    ) +
    ggplot2::labs(title = title, x = xlab, y = NULL, color = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = plot.title.size, face = "bold"),
      axis.title.x = ggplot2::element_text(size = axis.title.x.size),
      axis.text.x  = ggplot2::element_text(size = axis.text.x.size),
      axis.text.y  = ggplot2::element_text(size = axis.text.y.size, face = "italic"),
      plot.margin  = ggplot2::margin(t = margin.t, r = margin.r,
                                     b = margin.b, l = left.margin),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.box = legend.direction,
      legend.box.just = "center",
      legend.text = ggplot2::element_text(size = legend.text.size),
      legend.title = ggplot2::element_text(size = legend.title.size)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        direction = legend.direction,
        barwidth  = legend.barwidth,
        barheight = legend.barheight,
        ticks = FALSE
      )
    )
}
