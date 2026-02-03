# =============================================================================
# Built-in Design Templates
# =============================================================================
# Pre-configured designs for common CRT types. Users can call these directly
# without writing any model or cost functions.

library(glmmrBase)
library(Matrix)

# =============================================================================
# Parallel CRT Design
# =============================================================================

#' Create a parallel cluster randomised trial design
#'
#' Two-arm parallel design where clusters are randomised to treatment or control.
#' Stage 1 collects baseline data, Stage 2 can add more observations and/or clusters.
#'
#' @param icc Intra-cluster correlation coefficient
#' @param cac Cluster autocorrelation (correlation between periods within cluster)
#' @param delta Treatment effect on linear predictor scale
#' @param k1 Stage 1 clusters per arm. Can be scalar or vector for exploration.
#' @param m1 Stage 1 individuals per cluster-period. Can be scalar or vector.
#' @param k2 Stage 2 new clusters per arm. Vector of values to search over.
#' @param m2 Stage 2 individuals per cluster-period. Vector or function of m1.
#' @param rho Ratio of cluster cost to individual cost (default 30)
#' @param n_arms Number of arms (default 2)
#'
#' @return A crt_design object
#'
#' @examples
#' design <- parallel_crt(
#'   icc = 0.05, cac = 0.8, delta = 0.25,
#'   k1 = c(8, 10, 12), m1 = seq(10, 50, by = 10),
#'   k2 = 0:4, rho = 30
#' )
#' results <- adaptive_analysis(design, target_power = 0.8)
#'
#' @export
parallel_crt <- function(icc, cac = 0.8, delta,
                         k1, m1,
                         k2 = 0:4, m2 = NULL,
                         rho = 30, n_arms = 2) {

  # Create specification
  spec <- crt_design_spec(
    stage1_params = c("k1", "m1"),
    stage2_params = c("k2", "m2"),
    resources = list(
      n_s1 = ~ n_arms * k1 * m1,
      n_s2 = ~ n_arms * k1 * m2 + n_arms * k2 * m2,
      clusters_s1 = ~ n_arms * k1,
      clusters_s2 = ~ n_arms * k2
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho"),
      stage2_resources = c("n_s2", "clusters_s2")
    ),
    model_builder = parallel_model_builder,
    n_arms = n_arms,
    design_type = "parallel"
  )

  # Default m2 range based on m1

if (is.null(m2)) {
    m2_fn <- function(s1) {
      m1_val <- s1$m1
      seq(max(5, floor(0.5 * m1_val)), ceiling(2 * m1_val), by = 5)
    }
  } else if (is.function(m2)) {
    m2_fn <- m2
  } else {
    m2_fn <- function(s1) m2
  }

  # Stage 2 grid function
  stage2_grid_fn <- function(s1) {
    expand.grid(m2 = m2_fn(s1), k2 = k2)
  }

  structure(
    list(
      spec = spec,
      fixed_params = list(icc = icc, cac = cac, delta = delta),
      stage1_grid = expand.grid(k1 = k1, m1 = m1),
      stage2_grid_fn = stage2_grid_fn,
      rho = rho
    ),
    class = c("parallel_crt", "crt_design")
  )
}

#' Model builder for parallel CRT (internal)
#' @keywords internal
parallel_model_builder <- function(design_params, fixed_params) {

  k1 <- fixed_params$k1
  m1 <- fixed_params$m1
  k2 <- design_params$k2
  m2 <- design_params$m2
  m2b <- design_params$m2b %||% design_params$m2
  icc <- fixed_params$icc
  cac <- fixed_params$cac
  delta <- fixed_params$delta
  cache <- fixed_params$cache

  # Stage 1 design
  df1 <- nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
  df1$int <- df1$int - 1
  df1$n <- m1
  df1$t <- 0

  # Stage 2: existing clusters
  df2 <- df1
  df2$n <- m2
  df2$t <- 1

  # Stage 2: new clusters (if any)
  if (k2 > 0) {
    df2b <- nelder(as.formula(paste0("~ int(2) > cl(", k2, ")")))
    df2b$int <- df2b$int - 1
    df2b$cl <- df2b$cl + k1 * 2
    df2b$n <- m2b
    df2b$t <- 1
    df <- rbind(df1, df2, df2b)
  } else {
    df <- rbind(df1, df2)
  }

  beta <- c(0, delta)
  v1 <- c(icc * cac, icc * (1 - cac))

  # Cache key
  cache_key <- paste(k1, k2, sep = "_")

  # Build or update model
  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- Model$new(
      ~ int + t + (1|gr(cl)) + (1|gr(cl,t)),
      data = df,
      family = gaussian(),
      mean = c(beta, 0),
      covariance = v1,
      weights = df$n,
      var_par = 1 - icc
    )
    cache$cache_key <- cache_key
  } else {
    glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    cache$mod$update_parameters(cov.pars = v1)
  }

  # Extract matrices
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  if (k2 > 0) n2 <- n2 + nrow(df2b)

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  idx1 <- 1:n1
  idx2 <- (n1 + 1):(n1 + n2)

  # Efficient score decomposition
  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = 2)

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = delta,
    n1 = n1,
    n2 = n2,
    m1 = m1,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2
  )
}

#' Fixed parallel model builder (for comparison)
#' @keywords internal
fixed_parallel_model_builder <- function(design_params, fixed_params) {

  k <- design_params$k
  m <- design_params$m
  icc <- fixed_params$icc
  cac <- fixed_params$cac %||% 0.8
  delta <- fixed_params$delta

  df <- nelder(as.formula(paste0("~ int(2) > cl(", k, ")")))
  df$int <- df$int - 1
  df$n <- m

  mod <- Model$new(
    ~ int + (1|gr(cl)),
    data = df,
    family = gaussian(),
    mean = c(0, delta),
    covariance = icc,
    weights = df$n,
    var_par = 1 - icc
  )

  X <- mod$mean$X
  D <- mod$covariance$D
  Z <- mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  V_inv_X <- solve(V, X)
  I_full <- crossprod(X, V_inv_X)
  I_eff <- 1 / solve(I_full)[2, 2]

  list(
    I_eff = I_eff,
    b1 = delta,
    k = k,
    m = m,
    n_total = 2 * k * m
  )
}

# =============================================================================
# Crossover CRT Design
# =============================================================================

#' Create a crossover cluster randomised trial design
#'
#' Two-period crossover design where clusters cross over between treatment
#' and control. Stage 2 can add more observations and/or clusters.
#'
#' @param icc Intra-cluster correlation coefficient
#' @param cac Cluster autocorrelation
#' @param delta Treatment effect on linear predictor scale
#' @param k1 Stage 1 clusters per arm
#' @param m1 Stage 1 individuals per cluster-period
#' @param k2 Stage 2 new clusters per arm (vector)
#' @param m2 Stage 2 individuals per cluster-period (vector or function)
#' @param rho Ratio of cluster cost to individual cost
#' @param n_arms Number of arms (default 2)
#'
#' @return A crt_design object
#' @export
crossover_crt <- function(icc, cac = 0.8, delta,
                          k1, m1,
                          k2 = 0:4, m2 = NULL,
                          rho = 30, n_arms = 2) {

  spec <- crt_design_spec(
    stage1_params = c("k1", "m1"),
    stage2_params = c("k2", "m2"),
    resources = list(
      n_s1 = ~ n_arms * k1 * m1 * 2,        # 2 periods in stage 1
      n_s2 = ~ n_arms * k1 * m2 + n_arms * k2 * m2 * 2,  # existing + new clusters
      clusters_s1 = ~ n_arms * k1,
      clusters_s2 = ~ n_arms * k2
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho"),
      stage2_resources = c("n_s2", "clusters_s2")
    ),
    model_builder = crossover_model_builder,
    n_arms = n_arms,
    design_type = "crossover"
  )

  if (is.null(m2)) {
    m2_fn <- function(s1) seq(max(5, floor(0.5 * s1$m1)), ceiling(2 * s1$m1), by = 5)
  } else if (is.function(m2)) {
    m2_fn <- m2
  } else {
    m2_fn <- function(s1) m2
  }

  stage2_grid_fn <- function(s1) {
    expand.grid(m2 = m2_fn(s1), k2 = k2)
  }

  structure(
    list(
      spec = spec,
      fixed_params = list(icc = icc, cac = cac, delta = delta),
      stage1_grid = expand.grid(k1 = k1, m1 = m1),
      stage2_grid_fn = stage2_grid_fn,
      rho = rho
    ),
    class = c("crossover_crt", "crt_design")
  )
}

#' Model builder for crossover CRT
#' @keywords internal
crossover_model_builder <- function(design_params, fixed_params) {

  k1 <- fixed_params$k1
  m1 <- fixed_params$m1
  k2 <- design_params$k2
  m2 <- design_params$m2
  m2b <- design_params$m2b %||% design_params$m2
  icc <- fixed_params$icc
  cac <- fixed_params$cac
  delta <- fixed_params$delta
  cache <- fixed_params$cache

  # Stage 1: crossover (2 periods)
  df1 <- nelder(as.formula(paste0("~ int(2) > cl(", k1, ") > t(2)")))
  df1$int <- df1$int - 1
  df1$t <- df1$t - 1
  # Treatment: arm 1 gets trt in period 2, arm 2 gets trt in period 1
  df1$trt <- ifelse((df1$int == 0 & df1$t == 1) | (df1$int == 1 & df1$t == 0), 1, 0)
  df1$n <- m1

  # Stage 2 period (period 3)
  df2_base <- nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
  df2_base$int <- df2_base$int - 1
  df2_base$t <- 2
  df2_base$trt <- df2_base$int  # Same as their original arm
  df2_base$n <- m2

  if (k2 > 0) {
    # New clusters get both periods
    df2_new <- nelder(as.formula(paste0("~ int(2) > cl(", k2, ") > t(2)")))
    df2_new$int <- df2_new$int - 1
    df2_new$cl <- df2_new$cl + k1 * 2
    df2_new$t <- df2_new$t + 1  # Periods 2 and 3
    df2_new$trt <- ifelse((df2_new$int == 0 & df2_new$t == 3) |
                            (df2_new$int == 1 & df2_new$t == 2), 1, 0)
    df2_new$n <- m2b
    df <- rbind(df1, df2_base, df2_new)
  } else {
    df <- rbind(df1, df2_base)
  }

  v1 <- c(icc * cac, icc * (1 - cac))
  cache_key <- paste(k1, k2, sep = "_")

  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- Model$new(
      ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
      data = df,
      family = gaussian(),
      mean = c(rep(0, length(unique(df$t))), delta, 0),
      covariance = v1,
      weights = df$n,
      var_par = 1 - icc
    )
    cache$cache_key <- cache_key
  } else {
    glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    cache$mod$update_parameters(cov.pars = v1)
  }

  n1 <- nrow(df1)
  n2 <- nrow(df) - n1

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  idx1 <- 1:n1
  idx2 <- (n1 + 1):nrow(df)

  # Treatment column is the last fixed effect before random effects
  j_trt <- ncol(X) - 1  # May need adjustment based on model structure

  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = delta,
    n1 = n1,
    n2 = n2,
    m1 = m1,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2
  )
}

# =============================================================================
# Stepped-Wedge / Rollout Design
# =============================================================================

#' Create a stepped-wedge (rollout) cluster randomised trial design
#'
#' Stepped-wedge design where clusters progressively switch from control to
#' treatment over time. Stage 2 can add more time periods with different
#' rollout rates.
#'
#' @param icc Intra-cluster correlation coefficient
#' @param cac Cluster autocorrelation
#' @param delta Treatment effect on linear predictor scale
#' @param K Total number of clusters (fixed throughout trial)
#' @param t1 Number of time periods in stage 1
#' @param r1 Clusters switching per period in stage 1
#' @param m Individuals per cluster-period (fixed)
#' @param t2 Stage 2 additional time periods (vector to search over)
#' @param r2 Stage 2 rollout rate - clusters per period (vector to search over)
#' @param rho Ratio of time period cost to individual cost
#'
#' @return A crt_design object
#' @export
stepped_wedge_crt <- function(icc, cac = 0.8, delta,
                               K, t1, r1 = NULL, m,
                               t2 = 1:6, r2 = NULL,
                               rho = 100) {

  # Default r1: spread evenly across stage 1
  if (is.null(r1)) r1 <- ceiling(K / t1)

  # Default r2: same options as reasonable spread
  if (is.null(r2)) r2 <- unique(c(1, 2, ceiling(K / 6), ceiling(K / 4), ceiling(K / 2)))
  r2 <- r2[r2 > 0 & r2 <= K]

  spec <- crt_design_spec(
    stage1_params = c("K", "t1", "r1", "m"),
    stage2_params = c("t2", "r2"),
    resources = list(
      n_s1 = ~ K * t1 * m,
      n_s2 = ~ K * t2 * m,
      periods_s1 = ~ t1,
      periods_s2 = ~ t2
    ),
    cost_structure = list(
      weights = c(n = 1, periods = "rho"),
      stage2_resources = c("n_s2", "periods_s2")
    ),
    model_builder = stepped_wedge_model_builder,
    n_arms = 1,  # Not really arms in SW
    design_type = "stepped_wedge"
  )

  # For SW, stage 1 grid might just be different K or t1 values
  if (length(K) > 1 || length(t1) > 1) {
    stage1_grid <- expand.grid(K = K, t1 = t1, r1 = r1[1], m = m)
  } else {
    stage1_grid <- data.frame(K = K, t1 = t1, r1 = r1, m = m)
  }

  stage2_grid_fn <- function(s1) {
    expand.grid(t2 = t2, r2 = r2)
  }

  structure(
    list(
      spec = spec,
      fixed_params = list(icc = icc, cac = cac, delta = delta),
      stage1_grid = stage1_grid,
      stage2_grid_fn = stage2_grid_fn,
      rho = rho
    ),
    class = c("stepped_wedge_crt", "crt_design")
  )
}

#' Model builder for stepped-wedge design
#' @keywords internal
stepped_wedge_model_builder <- function(design_params, fixed_params) {

  K <- fixed_params$K
  t1 <- fixed_params$t1
  r1 <- fixed_params$r1
  m <- fixed_params$m
  t2 <- design_params$t2
  r2 <- design_params$r2
  icc <- fixed_params$icc
  cac <- fixed_params$cac
  delta <- fixed_params$delta
  cache <- fixed_params$cache

  T_total <- t1 + t2

  # Calculate switch times for each cluster
  switch_times <- numeric(K)
  clusters_assigned <- 0

  # Stage 1 assignment
  for (t in 1:t1) {
    n_switch <- min(r1, K - clusters_assigned)
    if (n_switch > 0) {
      idx <- (clusters_assigned + 1):(clusters_assigned + n_switch)
      switch_times[idx] <- t
      clusters_assigned <- clusters_assigned + n_switch
    }
  }

  # Stage 2 assignment
  if (clusters_assigned < K && t2 > 0) {
    for (t in (t1 + 1):T_total) {
      n_switch <- min(r2, K - clusters_assigned)
      if (n_switch > 0) {
        idx <- (clusters_assigned + 1):(clusters_assigned + n_switch)
        switch_times[idx] <- t
        clusters_assigned <- clusters_assigned + n_switch
      }
      if (clusters_assigned >= K) break
    }
  }

  # Remaining clusters switch at last period
  if (any(switch_times == 0)) {
    switch_times[switch_times == 0] <- T_total
  }

  # Build data frame
  df <- expand.grid(cl = 1:K, t = 1:T_total)
  df$trt <- as.integer(df$t >= switch_times[df$cl])
  df$n <- m

  idx1 <- which(df$t <= t1)
  idx2 <- which(df$t > t1)

  v1 <- c(icc * cac, icc * (1 - cac))
  cache_key <- paste(K, T_total, r2, sep = "_")

  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- Model$new(
      as.formula("~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t))"),
      data = df,
      family = gaussian(),
      mean = c(rep(0, T_total), delta, 0),
      covariance = v1,
      weights = df$n,
      var_par = 1 - icc
    )
    cache$cache_key <- cache_key
  } else {
    glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    cache$mod$update_parameters(cov.pars = v1)
  }

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  # Treatment is the last fixed effect column
  j_trt <- ncol(X) - 1

  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)

  n1 <- length(idx1) * m
  n2 <- length(idx2) * m

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = delta,
    n1 = n1,
    n2 = n2,
    K = K,
    t1 = t1,
    t2 = t2,
    r1 = r1,
    r2 = r2,
    m = m,
    T_total = T_total
  )
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
