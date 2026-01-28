waterfall.svcf.cate <-
function(svcf, tau, X = NULL, i = 1, k = 10,
                                base.value = NULL,
                                title = NULL,
                                xlab = "Treatment Effect (Cumulative)",
                                pos.col = "#E07A00",
                                neg.col = "#10978F",
                                rect.alpha = 0.85,
                                seg.alpha = 0.80,
                                seg.size = 0.9,
                                rect.h = 0.35,
                                base.linetype = "dashed",
                                base.alpha = 0.5,
                                digits = 3,
                                value.digits = 5,
                                value.thresh = 1e-5,
                                plot.title.size = 15,
                                axis.title.x.size = 10,
                                axis.text.x.size = 10,
                                axis.text.y.size = 12,
                                left.margin = 8,
                                margin.t = 6, margin.r = 6, margin.b = 6) {
  
  if (is.null(base.value)) base.value <- mean(tau, na.rm = TRUE)
  
  stopifnot(is.matrix(svcf) || is.data.frame(svcf))
  stopifnot(!is.null(colnames(svcf)))
  stopifnot(i >= 1, i == as.integer(i), i <= nrow(svcf))
  stopifnot(k >= 1, k == as.integer(k))
  stopifnot(is.numeric(base.value), length(base.value) == 1, is.finite(base.value))
  stopifnot(digits >= 0, digits == as.integer(digits))
  stopifnot(value.digits >= 0, value.digits == as.integer(value.digits))
  stopifnot(rect.alpha >= 0, rect.alpha <= 1,
            seg.alpha  >= 0, seg.alpha  <= 1,
            base.alpha >= 0, base.alpha <= 1)
  
  if (!is.null(X)) {
    stopifnot(is.matrix(X) || is.data.frame(X))
    stopifnot(i <= nrow(X))
    if (is.null(colnames(X))) stop("If 'X' is provided, it must have column names.")
  }
  
  fmt.value <- function(z) {
    if (is.na(z)) return(NA_character_)
    z2 <- round(as.numeric(z), value.digits + 1)
    if (!is.finite(z2)) return(NA_character_)
    if (abs(z2) < value.thresh) {
      paste0("<", formatC(value.thresh, format = "f", digits = value.digits))
    } else {
      formatC(z2, format = "f", digits = value.digits)
    }
  }
  
  s <- as.numeric(svcf[i, ])
  names(s) <- colnames(svcf)
  
  ord  <- order(abs(s), decreasing = TRUE)
  keep <- ord[seq_len(min(k, length(ord)))]
  rest <- setdiff(seq_along(s), keep)
  
  df <- data.frame(
    feature = names(s)[keep],
    shap = s[keep],
    stringsAsFactors = FALSE
  )
  
  if (length(rest) > 0) {
    df <- rbind(
      df,
      data.frame(feature = "Other", shap = sum(s[rest]), stringsAsFactors = FALSE)
    )
  }
  
  if (!is.null(X)) {
    xrow <- X[i, , drop = FALSE]
    df$value <- NA_real_
    
    for (j in seq_len(nrow(df))) {
      f <- df$feature[j]
      if (f != "Other" && f %in% colnames(xrow)) {
        df$value[j] <- as.numeric(xrow[[f]])
      }
    }
    
    df$label <- ifelse(
      is.na(df$value),
      df$feature,
      paste0(df$feature, "\n(", vapply(df$value, fmt.value, character(1)), ")")
    )
  } else {
    df$label <- df$feature
  }
  
  ## ---- cumulative construction
  df$start <- base.value + c(0, head(cumsum(df$shap), -1))
  df$end   <- base.value + cumsum(df$shap)
  
  df$xmin <- pmin(df$start, df$end)
  df$xmax <- pmax(df$start, df$end)
  df$dir  <- ifelse(df$shap >= 0, "positive", "negative")
  
  df$label <- factor(df$label, levels = rev(df$label))
  
  ggplot2::ggplot(df, ggplot2::aes(y = label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = start, xend = end, yend = label),
      linewidth = seg.size, lineend = "round", alpha = seg.alpha
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = as.numeric(label) - rect.h,
        ymax = as.numeric(label) + rect.h,
        fill = dir
      ),
      alpha = rect.alpha
    ) +
    ggplot2::geom_vline(
      xintercept = base.value,
      linetype = base.linetype,
      alpha = base.alpha
    ) +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    ggplot2::scale_x_continuous(
      labels = function(x) formatC(x, format = "f", digits = digits)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = plot.title.size, face = "bold"),
      axis.title.x = ggplot2::element_text(size = axis.title.x.size),
      axis.text.x  = ggplot2::element_text(size = axis.text.x.size),
      axis.text.y  = ggplot2::element_text(size = axis.text.y.size, face = "italic"),
      plot.margin  = ggplot2::margin(t = margin.t, r = margin.r,
                                     b = margin.b, l = left.margin),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::scale_fill_manual(values = c(positive = pos.col, negative = neg.col))
}
