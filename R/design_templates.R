# =============================================================================
# Built-in Design Templates
# =============================================================================
# Pre-configured designs for common CRT types. Users can call these directly
# without writing any model or cost functions.
#' Convert parameters to GLM scale for different families
#'
#' @param family Character: "gaussian", "binomial", or "poisson"
#' @param baseline Baseline outcome (probability for binomial, rate for poisson, mean for gaussian)
#' @param delta Treatment effect on natural scale
#' @param icc Intra-cluster correlation on natural scale
#' @param cac Cluster autocorrelation (if applicable)
#' @return List with beta, covariance, family object, and scale info
convert_to_glm_scale <- function(family = c("gaussian", "binomial", "poisson"),
                                 baseline = 0,
                                 delta,
                                 icc,
                                 cac = NULL) {


  family <- match.arg(family)

  if (family == "gaussian") {
    # Gaussian: parameters on natural scale
    beta <- c(baseline, delta)

    if (!is.null(cac)) {
      cov_pars <- c(icc * cac, icc * (1 - cac))
    } else {
      cov_pars <- icc
    }

    var_par <- 1 - icc
    family_obj <- gaussian()

    return(list(
      beta = beta,
      cov_pars = cov_pars,
      var_par = var_par,
      family = family_obj,
      b1_original = delta,  # Treatment effect on original scale
      link = "identity"
    ))

  } else if (family == "binomial") {
    # Binomial: convert to logit scale

    # Convert baseline and treatment effect to logit scale
    b0 <- log(baseline / (1 - baseline))  # logit(baseline)
    p1 <- baseline + delta  # probability under treatment

    # Check validity
    if (p1 <= 0 || p1 >= 1) {
      stop(sprintf("Treatment probability %.3f is outside (0,1). Adjust baseline or delta.", p1))
    }

    b1 <- log(p1 / (1 - p1)) - b0  # logit(p1) - logit(p0)
    beta <- c(b0, b1)

    # Convert ICC to variance on logit scale
    # Variance of Y on probability scale: p(1-p)
    # Variance of logit(Y): approximately 1 / (p(1-p)) for logistic
    v_prob <- baseline * (1 - baseline)  # variance on prob scale
    v_between_prob <- v_prob * (icc / (1 - icc))  # between-cluster variance on prob scale

    # Delta method: Var(logit(p)) ≈ Var(p) / (p(1-p))^2
    deriv <- 1 / (baseline * (1 - baseline))  # derivative of logit at baseline
    tau_sq <- v_between_prob * deriv^2

    if (!is.null(cac)) {
      cov_pars <- c(tau_sq * cac, tau_sq * (1 - cac))
    } else {
      cov_pars <- tau_sq
    }

    var_par <- 1  # For binomial, residual variance is determined by p(1-p)/n
    family_obj <- binomial()

    return(list(
      beta = beta,
      cov_pars = cov_pars,
      var_par = var_par,
      family = family_obj,
      b1_original = delta,  # Treatment effect on probability scale
      b1_link = b1,         # Treatment effect on logit scale (log-OR)
      link = "logit"
    ))

  } else if (family == "poisson") {
    # Poisson: convert to log scale

    # baseline is the rate under control
    b0 <- log(baseline)
    rate1 <- baseline + delta  # rate under treatment

    # Check validity
    if (rate1 <= 0) {
      stop(sprintf("Treatment rate %.3f must be positive. Adjust baseline or delta.", rate1))
    }

    b1 <- log(rate1) - b0  # log(rate1) - log(rate0) = log rate ratio
    beta <- c(b0, b1)

    # Convert ICC to variance on log scale
    # For Poisson, variance = mean, so ICC relates to overdispersion
    # Var(log(Y)) ≈ 1/mu for Poisson
    # Between-cluster variance on log scale
    v_log <- 1 / baseline  # approximate variance of log(Y)
    v_between_log <- v_log * (icc / (1 - icc))
    tau_sq <- v_between_log

    if (!is.null(cac)) {
      cov_pars <- c(tau_sq * cac, tau_sq * (1 - cac))
    } else {
      cov_pars <- tau_sq
    }

    var_par <- 1  # For Poisson, residual is determined by mean
    family_obj <- poisson()

    return(list(
      beta = beta,
      cov_pars = cov_pars,
      var_par = var_par,
      family = family_obj,
      b1_original = delta,  # Treatment effect on rate scale
      b1_link = b1,         # Treatment effect on log scale (log-RR)
      link = "log"
    ))
  }
}
# =============================================================================
# Parallel CRT Design
# =============================================================================

#' Create a parallel cluster randomised trial design
#'
#' Two-arm parallel design where clusters are randomised to treatment or control.
#' Stage 2 can add more observations and/or clusters.
#'
#' @param icc Intra-cluster correlation coefficient
#' @param delta Treatment effect (on natural scale)
#' @param k1 Stage 1 clusters per arm (scalar or vector)
#' @param m1 Stage 1 individuals per cluster (scalar or vector)
#' @param k2 Stage 2 new clusters per arm (vector, default 0:4)
#' @param m2 Stage 2 individuals per cluster (vector or function, default based on m1)
#' @param rho Ratio of cluster cost to individual cost
#' @param cac Cluster autocorrelation (default 0.8)
#' @param family Distribution family: "gaussian", "binomial", or "poisson"
#' @param baseline Baseline outcome (required for binomial/poisson)
#'
#' @return A crt_design object
#' @export
parallel_crt <- function(icc, delta, k1, m1,
                         k2 = 0:4, m2 = NULL,
                         rho = 1, cac = 0.8,
                         family = "gaussian", baseline = NULL) {

  # Validate family
  family <- match.arg(family, c("gaussian", "binomial", "poisson"))

  # Baseline required for non-Gaussian
  if (family != "gaussian" && is.null(baseline)) {
    stop("baseline is required for ", family, " family")
  }

  # Validate baseline for binomial
  if (family == "binomial") {
    if (baseline <= 0 || baseline >= 1) {
      stop("baseline must be between 0 and 1 for binomial family")
    }
    if (baseline + delta <= 0 || baseline + delta >= 1) {
      stop("baseline + delta must be between 0 and 1 for binomial family")
    }
  }

  # Validate baseline for poisson
  if (family == "poisson") {
    if (baseline <= 0) {
      stop("baseline must be positive for poisson family")
    }
    if (baseline + delta <= 0) {
      stop("baseline + delta must be positive for poisson family")
    }
  }

  n_arms <- 2

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

  # Stage 2 grid function
  if (is.null(m2)) {
    m2_fn <- function(s1) seq(max(5, floor(0.5 * s1$m1)), ceiling(2 * s1$m1), by = 10)
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
      fixed_params = list(
        icc = icc,
        cac = cac,
        delta = delta,
        family = family,
        baseline = baseline %||% 0
      ),
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

  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0

  # Convert parameters to appropriate scale
  glm_params <- convert_to_glm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = cac
  )

  # Stage 1 design
  df1 <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
  df1$int <- df1$int - 1
  df1$n <- m1
  df1$t <- 0

  # Stage 2: existing clusters
  df2 <- df1
  df2$n <- m2
  df2$t <- 1

  # Stage 2: new clusters (if any)
  if (k2 > 0) {
    df2b <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k2, ")")))
    df2b$int <- df2b$int - 1
    df2b$cl <- df2b$cl + k1 * 2
    df2b$n <- m2b
    df2b$t <- 1
    df <- rbind(df1, df2, df2b)
  } else {
    df <- rbind(df1, df2)
  }

  # Cache key includes family
  cache_key <- paste(k1, k2, family, sep = "_")

  # Build or update model
  if (is.null(cache$mod) || cache$cache_key != cache_key) {

    # Build model based on family
    if (family == "gaussian") {
      cache$mod <- glmmrBase::Model$new(
        ~ int + t + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, 0),
        covariance = glm_params$cov_pars,
        weights = df$n,
        var_par = glm_params$var_par
      )
    } else if (family == "binomial") {
      cache$mod <- glmmrBase::Model$new(
        ~ int + t + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, 0),
        covariance = glm_params$cov_pars,
        trials = df$n
      )
    } else if (family == "poisson") {
      # For Poisson, n represents exposure (e.g., person-time)
      cache$mod <- glmmrBase::Model$new(
        ~ int + t + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, 0),
        covariance = glm_params$cov_pars,
        offset = log(df$n)
      )
    }

    cache$cache_key <- cache_key

  } else {
    # Update existing model
    if (family == "gaussian") {
      glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    } else if (family == "binomial") {
      glmmrBase:::Model__set_trials(cache$mod$.__enclos_env__$private$ptr, as.integer(df$n))
    } else if (family == "poisson") {
      glmmrBase:::Model__set_offset(cache$mod$.__enclos_env__$private$ptr, log(df$n))
    }
    cache$mod$update_parameters(cov.pars = glm_params$cov_pars)
  }

  # Extract matrices
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  if (k2 > 0) n2 <- n2 + nrow(df2b)

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)

  # For non-Gaussian, need to account for variance function
  S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  idx1 <- 1:n1
  idx2 <- (n1 + 1):(n1 + n2)

  # Efficient score decomposition
  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = 2)

  # Degrees of freedom
  n_cluster_periods <- 2 * (k1 + k2)
  n_fixed <- 2  # intercept + treatment
  df_full <- n_cluster_periods - n_fixed
  n_cluster_periods_s1 <- 2 * k1
  df_s1 <- n_cluster_periods_s1 - n_fixed

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = glm_params$b1_link %||% glm_params$beta[2],  # Effect on link scale
    b1_original = glm_params$b1_original,              # Effect on original scale
    n1 = n1,
    n2 = n2,
    m1 = m1,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2,
    df_s1 = df_s1,
    df_full = df_full,
    family = family,
    link = glm_params$link
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

  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0

  # Convert parameters
  glm_params <- convert_to_glm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = NULL  # No CAC for simple fixed design
  )

  # Two-arm parallel design
  n_arms <- 2

  # Build data
  df <- data.frame(
    cl = rep(1:(n_arms * k), each = 1),
    trt = rep(0:1, each = k)
  )
  df$n <- m

  # Build model based on family
  if (family == "gaussian") {
    mod <- glmmrBase::Model$new(
      formula = ~ trt + (1|gr(cl)),
      data = df,
      family = glm_params$family,
      mean = glm_params$beta,
      covariance = glm_params$cov_pars,
      weights = df$n,
      var_par = glm_params$var_par
    )
  } else if (family == "binomial") {
    mod <- glmmrBase::Model$new(
      formula = ~ trt + (1|gr(cl)),
      data = df,
      family = glm_params$family,
      mean = glm_params$beta,
      covariance = glm_params$cov_pars,
      trials = df$n
    )
  } else if (family == "poisson") {
    mod <- glmmrBase::Model$new(
      formula = ~ trt + (1|gr(cl)),
      data = df,
      family = glm_params$family,
      mean = glm_params$beta,
      covariance = glm_params$cov_pars,
      offset = log(df$n)
    )
  }

  X <- mod$mean$X

  # Get covariance matrix
  D <- mod$covariance$D
  Z <- mod$covariance$Z
  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  # Full information
  V_inv_X <- solve(V, X)
  I_full <- crossprod(X, V_inv_X)

  # Treatment is column 2
  j_trt <- 2
  I_eff <- 1 / solve(I_full)[j_trt, j_trt]

  # Total sample size
  n_total <- n_arms * k * m

  list(
    I_eff = I_eff,
    b1 = glm_params$b1_link %||% glm_params$beta[2],
    b1_original = glm_params$b1_original,
    k = k,
    m = m,
    n_total = n_total,
    family = family,
    link = glm_params$link
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
#' @param delta Treatment effect (on natural scale)
#' @param k1 Stage 1 clusters per arm
#' @param m1 Stage 1 individuals per cluster-period
#' @param k2 Stage 2 new clusters per arm (vector)
#' @param m2 Stage 2 individuals per cluster-period (vector or function)
#' @param rho Ratio of cluster cost to individual cost
#' @param n_arms Number of arms (default 2)
#' @param family Distribution family: "gaussian", "binomial", or "poisson"
#' @param baseline Baseline outcome (required for binomial/poisson)
#'
#' @return A crt_design object
#' @export
crossover_crt <- function(icc, cac = 0.8, delta,
                          k1, m1,
                          k2 = 0:4, m2 = NULL,
                          rho = 30, n_arms = 2,
                          family = "gaussian", baseline = NULL) {

  # Validate family
  family <- match.arg(family, c("gaussian", "binomial", "poisson"))

  # Baseline required for non-Gaussian
  if (family != "gaussian" && is.null(baseline)) {
    stop("baseline is required for ", family, " family")
  }

  # Validate baseline for binomial
  if (family == "binomial") {
    if (baseline <= 0 || baseline >= 1) {
      stop("baseline must be between 0 and 1 for binomial family")
    }
    if (baseline + delta <= 0 || baseline + delta >= 1) {
      stop("baseline + delta must be between 0 and 1 for binomial family")
    }
  }

  # Validate baseline for poisson
  if (family == "poisson") {
    if (baseline <= 0) {
      stop("baseline must be positive for poisson family")
    }
    if (baseline + delta <= 0) {
      stop("baseline + delta must be positive for poisson family")
    }
  }

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
      fixed_params = list(
        icc = icc,
        cac = cac,
        delta = delta,
        family = family,
        baseline = baseline %||% 0
      ),
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

  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0

  # Convert parameters to appropriate scale
  glm_params <- convert_to_glm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = cac
  )

  # Stage 1: crossover (2 periods)
  df1 <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ") > t(2)")))
  df1$int <- df1$int - 1
  df1$t <- df1$t - 1
  # Treatment: arm 1 gets trt in period 2, arm 2 gets trt in period 1
  df1$trt <- ifelse((df1$int == 0 & df1$t == 1) | (df1$int == 1 & df1$t == 0), 1, 0)
  df1$n <- m1

  # Stage 2 period (period 3)
  df2_base <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
  df2_base$int <- df2_base$int - 1
  df2_base$t <- 2
  df2_base$trt <- df2_base$int  # Same as their original arm
  df2_base$n <- m2

  if (k2 > 0) {
    # New clusters get both periods
    df2_new <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k2, ") > t(2)")))
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

  # Number of time periods
  n_periods <- length(unique(df$t))

  # Mean parameters: intercept + (n_periods - 1) time effects + treatment + extra 0 for random
  # glm_params$beta = c(intercept, treatment)
  # Need: c(intercept, time2, time3, ..., treatment, 0)
  beta_full <- c(glm_params$beta[1],           # Intercept
                 rep(0, n_periods - 1),         # Time effects
                 glm_params$beta[2]         )   # Treatment

  # Cache key includes family
  cache_key <- paste(k1, k2, family, sep = "_")

  if (is.null(cache$mod) || cache$cache_key != cache_key) {

    if (family == "gaussian") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = glm_params$cov_pars,
        weights = df$n,
        var_par = glm_params$var_par
      )
    } else if (family == "binomial") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = glm_params$cov_pars,
        trials = df$n
      )
    } else if (family == "poisson") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = glm_params$cov_pars,
        offset = log(df$n)
      )
    }

    cache$cache_key <- cache_key

  } else {
    # Update existing model
    if (family == "gaussian") {
      glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    } else if (family == "binomial") {
      glmmrBase:::Model__set_trials(cache$mod$.__enclos_env__$private$ptr, as.integer(df$n))
    } else if (family == "poisson") {
      glmmrBase:::Model__set_offset(cache$mod$.__enclos_env__$private$ptr, log(df$n))
    }
    cache$mod$update_parameters(cov.pars = glm_params$cov_pars)
  }

  n1 <- nrow(df1)
  n2 <- nrow(df) - n1

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  idx1 <- 1:n1
  idx2 <- (n1 + 1):nrow(df)

  # Treatment column - find it by name
  j_trt <- which(colnames(X) == "trt")
  if (length(j_trt) == 0) {
    # Fallback: last column before random effects
    j_trt <- ncol(X) - 1
  }

  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)

  # Degrees of freedom
  # Cluster-periods: stage 1 has 2*k1 clusters * 2 periods, stage 2 adds k1 + 2*k2 cluster-periods
  n_cluster_periods_s1 <- 2 * k1 * 2  # 2 arms * k1 clusters * 2 periods
  n_cluster_periods_full <- n_cluster_periods_s1 + 2 * k1 + 2 * k2 * 2  # + stage 2

  # Fixed effects: intercept + (n_periods - 1) time effects + treatment
  n_fixed <- 1 + (n_periods - 1) + 1

  df_s1 <- n_cluster_periods_s1 - (1 + 1 + 1)  # Only 2 periods in stage 1, so fewer time effects

  df_full <- n_cluster_periods_full - n_fixed

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = glm_params$b1_link %||% glm_params$beta[2],
    b1_original = glm_params$b1_original,
    n1 = n1,
    n2 = n2,
    m1 = m1,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2,
    df_s1 = df_s1,
    df_full = df_full,
    family = family,
    link = glm_params$link
  )
}

#' Fixed crossover model builder (for comparison)
#' @keywords internal
fixed_crossover_model_builder <- function(design_params, fixed_params) {

  k <- design_params$k
  m <- design_params$m
  icc <- fixed_params$icc
  cac <- fixed_params$cac %||% 0.8
  delta <- fixed_params$delta

  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0

  # Convert parameters
  glm_params <- convert_to_glm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = cac
  )

  # Two-arm crossover design (2 periods)
  n_arms <- 2
  n_periods <- 2

  # Build data: each cluster has 2 periods
  df <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k, ") > t(2)")))
  df$int <- df$int - 1
  df$t <- df$t - 1
  # Treatment: arm 0 gets trt in period 1, arm 1 gets trt in period 0
  df$trt <- ifelse((df$int == 0 & df$t == 1) | (df$int == 1 & df$t == 0), 1, 0)
  df$n <- m

  # Mean parameters: intercept + time effect + treatment
  beta_full <- c(glm_params$beta[1],   # Intercept
                 0,                      # Time effect (period 2 vs 1)
                 glm_params$beta[2])     # Treatment

  # Build model based on family
  if (family == "gaussian") {
    mod <- glmmrBase::Model$new(
      formula = ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
      data = df,
      family = glm_params$family,
      mean = c(beta_full, 0),
      covariance = glm_params$cov_pars,
      weights = df$n,
      var_par = glm_params$var_par
    )
  } else if (family == "binomial") {
    mod <- glmmrBase::Model$new(
      formula = ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
      data = df,
      family = glm_params$family,
      mean = c(beta_full, 0),
      covariance = glm_params$cov_pars,
      trials = df$n
    )
  } else if (family == "poisson") {
    mod <- glmmrBase::Model$new(
      formula = ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
      data = df,
      family = glm_params$family,
      mean = c(beta_full, 0),
      covariance = glm_params$cov_pars,
      offset = log(df$n)
    )
  }

  X <- mod$mean$X
  D <- mod$covariance$D
  Z <- mod$covariance$Z

  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  # Full information
  V_inv_X <- solve(V, X)
  I_full <- crossprod(X, V_inv_X)

  # Treatment column
  j_trt <- which(colnames(X) == "trt")
  if (length(j_trt) == 0) {
    j_trt <- ncol(X) - 1
  }

  I_eff <- 1 / solve(I_full)[j_trt, j_trt]

  # Total sample size: clusters * periods * individuals
  n_total <- n_arms * k * n_periods * m

  list(
    I_eff = I_eff,
    b1 = glm_params$b1_link %||% glm_params$beta[2],
    b1_original = glm_params$b1_original,
    k = k,
    m = m,
    n_total = n_total,
    n_clusters = n_arms * k,
    n_periods = n_periods,
    family = family,
    link = glm_params$link
  )
}


# Null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
