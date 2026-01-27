svcf.sficf <-
function(tau, X, num.trees = 3000, num.trees.tune = 3000,
                       num.rep = 1, tune.oob = TRUE, mtry = NULL, min.node.size = NULL, seed = NULL, ...) {
  
  stopifnot(num.rep >= 1, num.rep == as.integer(num.rep))
  stopifnot(length(tau) == nrow(X))
  
  has_seed_arg <- "seed" %in% names(formals(ranger::ranger))
  
  if (num.rep == 1) {
    
    if (tune.oob) {
      dat <- data.frame(tau = tau, X)
      
      if (is.null(mtry)) {
        p <- ncol(X)
        mtry.grid <- unique(ceiling(c(sqrt(p), p/10, p/5, p/3)))
      } else {
        mtry.grid <- mtry
      }
      
      if (is.null(min.node.size)) {
        node.grid <- c(5L, 10L)
      } else {
        node.grid <- min.node.size
      }
      
      p <- ncol(X)
      mtry.grid <- as.integer(pmax(1L, pmin(p, mtry.grid)))
      node.grid <- as.integer(node.grid)
      
      grid <- expand.grid(mtry = mtry.grid, min.node.size = node.grid)
      oob.err <- numeric(nrow(grid))
      
      for (i in seq_len(nrow(grid))) {
        m  <- as.integer(grid$mtry[i])
        ns <- as.integer(grid$min.node.size[i])
        
        this.seed <- if (is.null(seed)) NULL else (seed + i - 1)
        
        if (!is.null(this.seed) && !has_seed_arg) {
          set.seed(this.seed)
        }
        
        if (is.null(this.seed) || !has_seed_arg) {
          rf.tmp <- ranger(
            tau ~ ., data = dat,
            mtry = m, min.node.size = ns,
            oob.error = TRUE,
            num.trees = num.trees.tune
          )
        } else {
          rf.tmp <- ranger(
            tau ~ ., data = dat,
            mtry = m, min.node.size = ns,
            oob.error = TRUE,
            num.trees = num.trees.tune,
            seed = this.seed
          )
        }
        
        oob.err[i] <- rf.tmp$prediction.error
      }
      
      res <- cbind(grid, oob.err = oob.err)
      
      idx <- which.min(res[, "oob.err"])
      best.m  <- as.integer(res[idx, "mtry"])
      best.ns <- as.integer(res[idx, "min.node.size"])
      
      this.seed <- seed
      
      if (!is.null(this.seed) && !has_seed_arg) {
        set.seed(this.seed)
      }
      
      if (is.null(this.seed) || !has_seed_arg) {
        rf <- ranger(
          tau ~ ., data = dat,
          mtry = best.m, min.node.size = best.ns,
          oob.error = FALSE, num.trees = num.trees,
          ...
        )
      } else {
        rf <- ranger(
          tau ~ ., data = dat,
          mtry = best.m, min.node.size = best.ns,
          oob.error = FALSE, num.trees = num.trees,
          seed = this.seed,
          ...
        )
      }
      
      uni <- ranger.unify(rf, X)
      ts <- treeshap(uni, X)
      
      svcf <- ts$shaps
      
      s_all <- rowSums(svcf)
      ove.var <- mean(s_all^2)
      ind.var <- numeric(ncol(svcf))
      
      for (i in seq_len(ncol(svcf))) {
        s_minus_i <- s_all - svcf[, i]
        ind.var[i] <- mean(s_minus_i^2)
      }
      
      sficf <- 1 - ind.var/ove.var
      sficf <- sficf/sum(sficf)
      names(sficf) <- colnames(X)
      
      out <- list(svcf = svcf, sficf = sficf, tau = tau, res = res)
      return(out)
    }
    
    if (!tune.oob) {
      dat <- data.frame(tau = tau, X)
      
      this.seed <- seed
      
      if (!is.null(this.seed) && !has_seed_arg) {
        set.seed(this.seed)
      }
      
      if (is.null(this.seed) || !has_seed_arg) {
        rf <- ranger(
          tau ~ ., data = dat,
          mtry = mtry,
          min.node.size = min.node.size,
          oob.error = FALSE, num.trees = num.trees,
          ...
        )
      } else {
        rf <- ranger(
          tau ~ ., data = dat,
          mtry = mtry,
          min.node.size = min.node.size,
          oob.error = FALSE, num.trees = num.trees,
          seed = this.seed,
          ...
        )
      }
      
      uni <- ranger.unify(rf, X)
      ts <- treeshap(uni, X)
      
      svcf <- ts$shaps
      
      s_all <- rowSums(svcf)
      ove.var <- mean(s_all^2)
      ind.var <- numeric(ncol(svcf))
      
      for (i in seq_len(ncol(svcf))) {
        s_minus_i <- s_all - svcf[, i]
        ind.var[i] <- mean(s_minus_i^2)
      }
      
      sficf <- 1 - ind.var/ove.var
      sficf <- sficf/sum(sficf)
      names(sficf) <- colnames(X)
      
      out <- list(svcf = svcf, sficf = sficf, tau = tau, res = NULL)
      return(out)
    }
  }
  
  if (num.rep > 1) {
    
    dat <- data.frame(tau = tau, X)
    
    if (tune.oob) {
      if (is.null(mtry)) {
        p <- ncol(X)
        mtry.grid <- unique(ceiling(c(sqrt(p), p/10, p/5, p/3)))
      } else {
        mtry.grid <- mtry
      }
      
      if (is.null(min.node.size)) {
        node.grid <- c(5L, 10L)
      } else {
        node.grid <- min.node.size
      }
      
      p <- ncol(X)
      mtry.grid <- as.integer(pmax(1L, pmin(p, mtry.grid)))
      node.grid <- as.integer(node.grid)
      
      grid <- expand.grid(mtry = mtry.grid, min.node.size = node.grid)
      oob.err <- numeric(nrow(grid))
      
      for (i in seq_len(nrow(grid))) {
        m  <- as.integer(grid$mtry[i])
        ns <- as.integer(grid$min.node.size[i])
        
        this.seed <- if (is.null(seed)) NULL else (seed + i - 1)
        
        if (!is.null(this.seed) && !has_seed_arg) {
          set.seed(this.seed)
        }
        
        if (is.null(this.seed) || !has_seed_arg) {
          rf.tmp <- ranger(
            tau ~ ., data = dat,
            mtry = m, min.node.size = ns,
            oob.error = TRUE, num.trees = num.trees.tune
          )
        } else {
          rf.tmp <- ranger(
            tau ~ ., data = dat,
            mtry = m, min.node.size = ns,
            oob.error = TRUE, num.trees = num.trees.tune,
            seed = this.seed
          )
        }
        
        oob.err[i] <- rf.tmp$prediction.error
      }
      
      res <- cbind(grid, oob.err = oob.err)
      idx <- which.min(res[, "oob.err"])
      best.m  <- as.integer(res[idx, "mtry"])
      best.ns <- as.integer(res[idx, "min.node.size"])
      
    } else {
      res <- NULL
    }
    
    n <- nrow(X)
    p <- ncol(X)
    
    svcf.sum  <- matrix(0, nrow = n, ncol = p)
    sficf.sum <- numeric(p)
    
    for (r in seq_len(num.rep)) {
      
      this.seed <- if (is.null(seed)) NULL else (seed + r - 1)
      
      if (!is.null(this.seed) && !has_seed_arg) {
        set.seed(this.seed)
      }
      
      if (tune.oob) {
        if (is.null(this.seed) || !has_seed_arg) {
          rf <- ranger(
            tau ~ ., data = dat,
            mtry = best.m, min.node.size = best.ns,
            oob.error = FALSE, num.trees = num.trees,
            ...
          )
        } else {
          rf <- ranger(
            tau ~ ., data = dat,
            mtry = best.m, min.node.size = best.ns,
            oob.error = FALSE, num.trees = num.trees,
            seed = this.seed,
            ...
          )
        }
      } else {
        if (is.null(this.seed) || !has_seed_arg) {
          rf <- ranger(
            tau ~ ., data = dat,
            mtry = mtry,
            min.node.size = min.node.size,
            oob.error = FALSE, num.trees = num.trees,
            ...
          )
        } else {
          rf <- ranger(
            tau ~ ., data = dat,
            mtry = mtry,
            min.node.size = min.node.size,
            oob.error = FALSE, num.trees = num.trees,
            seed = this.seed,
            ...
          )
        }
      }
      
      uni <- ranger.unify(rf, X)
      ts <- treeshap(uni, X)
      svcf.r <- ts$shaps
      
      s_all <- rowSums(svcf.r)
      ove.var <- mean(s_all^2)
      ind.var <- numeric(p)
      
      for (j in seq_len(p)) {
        s_minus_j <- s_all - svcf.r[, j]
        ind.var[j] <- mean(s_minus_j^2)
      }
      
      sficf.r <- 1 - ind.var / ove.var
      sficf.r <- sficf.r / sum(sficf.r)
      
      svcf.sum  <- svcf.sum + svcf.r
      sficf.sum <- sficf.sum + sficf.r
    }
    
    svcf  <- svcf.sum / num.rep
    sficf <- sficf.sum / num.rep
    colnames(svcf) <- colnames(X)
    names(sficf) <- colnames(X)
    
    out <- list(svcf = svcf, sficf = sficf, tau = tau, res = res)
    return(out)
  }
}
