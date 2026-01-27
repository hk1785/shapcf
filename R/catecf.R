catecf <-
function(Y, X, W,
                   num.trees = 3000, num.rep = 1, seed = NULL,
                   tune.parameters = "all",
                   compute.oob.predictions = TRUE, ...) {
  
  stopifnot(length(Y) == nrow(X), length(W) == length(Y))
  stopifnot(num.rep >= 1, num.rep == as.integer(num.rep))
  
  has_seed_arg <- "seed" %in% names(formals(grf::causal_forest))
  
  if (num.rep == 1) {
    
    this.seed <- seed
    
    if (!is.null(this.seed) && !has_seed_arg) {
      set.seed(this.seed)
    }
    
    if (is.null(this.seed) || !has_seed_arg) {
      fit <- causal_forest(
        Y = Y, X = X, W = W,
        tune.parameters = tune.parameters,
        compute.oob.predictions = compute.oob.predictions,
        num.trees = num.trees, ...
      )
    } else {
      fit <- causal_forest(
        Y = Y, X = X, W = W,
        tune.parameters = tune.parameters,
        compute.oob.predictions = compute.oob.predictions,
        num.trees = num.trees,
        seed = this.seed, ...
      )
    }
    
    out <- as.numeric(predict(fit)$predictions)
    return(out)
  }
  
  if (num.rep > 1) {
    
    preds <- matrix(NA_real_, nrow = length(Y), ncol = num.rep)
    
    for (r in seq_len(num.rep)) {
      
      this.seed <- if (is.null(seed)) NULL else (seed + r - 1)
      
      if (!is.null(this.seed) && !has_seed_arg) {
        set.seed(this.seed)
      }
      
      if (is.null(this.seed) || !has_seed_arg) {
        fit <- causal_forest(
          Y = Y, X = X, W = W,
          tune.parameters = tune.parameters,
          compute.oob.predictions = compute.oob.predictions,
          num.trees = num.trees, ...
        )
      } else {
        fit <- causal_forest(
          Y = Y, X = X, W = W,
          tune.parameters = tune.parameters,
          compute.oob.predictions = compute.oob.predictions,
          num.trees = num.trees,
          seed = this.seed, ...
        )
      }
      
      preds[,r] <- as.numeric(predict(fit)$predictions)
    }
    
    out <- rowMeans(preds)
    return(out)
  }
}
