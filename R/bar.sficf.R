bar.sficf <-
function(sficf, k = 10,
                      title = NULL,
                      xlab = "Feature Importance",
                      bar.col = "#6B7280",
                      bar.alpha = 0.9,
                      bar.width = 0.7,
                      base.size = 12,
                      plot.title.size = 12,
                      axis.title.x.size = 11,
                      axis.text.x.size = 10,
                      axis.text.y.size = 11,
                      left.margin = 8,
                      margin.t = 6, margin.r = 6, margin.b = 6) {
  
  stopifnot(is.numeric(sficf), !is.null(names(sficf)))
  stopifnot(k >= 1, k == as.integer(k))
  stopifnot(bar.alpha >= 0, bar.alpha <= 1)
  
  sficf <- sficf[is.finite(sficf)]
  if (length(sficf) == 0) stop("'sficf' has no finite values.")
  
  s2 <- sort(sficf, decreasing = TRUE)
  k2 <- min(k, length(s2))
  
  df <- data.frame(
    feature = names(s2)[seq_len(k2)],
    importance = as.numeric(s2[seq_len(k2)]),
    stringsAsFactors = FALSE
  )
  df$feature <- factor(df$feature, levels = rev(df$feature))
  
  ggplot2::ggplot(df, ggplot2::aes(x = importance, y = feature)) +
    ggplot2::geom_col(fill = bar.col, alpha = bar.alpha, width = bar.width) +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    ggplot2::theme_bw(base_size = base.size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = axis.title.x.size),
      axis.text.x  = ggplot2::element_text(size = axis.text.x.size),
      axis.text.y  = ggplot2::element_text(size = axis.text.y.size, face = "italic"),
      plot.title   = ggplot2::element_text(size = plot.title.size, face = "bold", hjust = 0.5),
      plot.margin  = ggplot2::margin(t = margin.t, r = margin.r, b = margin.b, l = left.margin)
    )
}
