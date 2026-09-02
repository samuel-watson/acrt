# =============================================================================
# Adaptive Cluster Trial Framework v2
# Generalised to support arbitrary stage 2 design parameters
# =============================================================================

## Density of Z1 = Phi^{-1}(F_{t_nu}(T1)) when T1 ~ t_nu(ncp = mu1).
## Under H0 (mu1 = 0) this is exactly dnorm(z); it reduces to dnorm(z - mu1)
## as nu -> Inf. Log-scale throughout, so the quadrature tails stay finite.
z1_density <- function(z, mu1, df) {
  if (!is.finite(df)) return(dnorm(z - mu1))
  tq <- suppressWarnings(ifelse(z >= 0,
               qt(pnorm(z, lower.tail = FALSE, log.p = TRUE), df = df,
                  lower.tail = FALSE, log.p = TRUE),
               qt(pnorm(z, log.p = TRUE), df = df, log.p = TRUE)))
  exp(dnorm(z, log = TRUE) +
        suppressWarnings(dt(tq, df = df, ncp = mu1, log = TRUE)) -
        suppressWarnings(dt(tq, df = df, log = TRUE)))
}

## P(|Z1| > efficacy_boundary). Two-sided in both branches so that it agrees
## with prob_efficacy by construction.
stage1_power <- function(efficacy_boundary, mu1, df_s1) {
  if (!is.finite(df_s1))
    return(pnorm(-efficacy_boundary - mu1) +
             pnorm(efficacy_boundary - mu1, lower.tail = FALSE))
  bt <- qt(pnorm(efficacy_boundary), df = df_s1)
  pt(bt, df = df_s1, ncp = mu1, lower.tail = FALSE) + pt(-bt, df = df_s1, ncp = mu1)
}

#' Compute efficient score decomposition for treatment effect
#'
#' Projects out nuisance parameters first, then decomposes by stage
#'
#' @param X Full design matrix
#' @param V Full marginal covariance matrix
#' @param idx1 Row indices for stage 1
#' @param idx2 Row indices for stage 2
#' @param j Column index for treatment effect (default 2)
#' @return List with I1_eff, I2_eff, I_eff (all scalars)
#' @export
efficient_score_decomposition <- function(X, V, idx1, idx2, j = 2) {

  n <- nrow(X)
  p <- ncol(X)

  # Partition design matrix
  X_trt <- X[, j, drop = FALSE]           # Treatment column (n x 1)
  X_nuis <- X[, -j, drop = FALSE]         # Nuisance columns (n x (p-1))

  # Solve V^{-1} times various quantities via Cholesky
  R <- chol(V)  # V = R'R

  # V^{-1} X_nuis
  V_inv_X_nuis <- backsolve(R, backsolve(R, X_nuis, transpose = TRUE))

  # V^{-1} X_trt
  V_inv_X_trt <- backsolve(R, backsolve(R, X_trt, transpose = TRUE))

  # (X_nuis' V^{-1} X_nuis)^{-1}
  I_nuis <- crossprod(X_nuis, V_inv_X_nuis)
  I_nuis_inv <- solve(I_nuis)

  # X_nuis' V^{-1} X_trt
  cross_info <- crossprod(X_nuis, V_inv_X_trt)

  # Projected treatment column: X_trt_tilde = X_trt - X_nuis * (I_nuis^{-1} * cross_info)
  adjustment <- X_nuis %*% (I_nuis_inv %*% cross_info)
  X_trt_tilde <- X_trt - adjustment  # n x 1 vector

  # Now decompose this 1D problem by stage
  X1_tilde <- X_trt_tilde[idx1, , drop = FALSE]
  X2_tilde <- X_trt_tilde[idx2, , drop = FALSE]

  V11 <- V[idx1, idx1, drop = FALSE]
  V12 <- V[idx1, idx2, drop = FALSE]
  V22 <- V[idx2, idx2, drop = FALSE]

  # Stage 1: I1_eff = X1_tilde' V11^{-1} X1_tilde
  R11 <- chol(V11)
  V11_inv_X1 <- backsolve(R11, backsolve(R11, X1_tilde, transpose = TRUE))
  I1_eff <- as.numeric(crossprod(X1_tilde, V11_inv_X1))

  # Schur complement: S = V22 - V21 V11^{-1} V12
  V11_inv_V12 <- backsolve(R11, backsolve(R11, V12, transpose = TRUE))
  S <- V22 - crossprod(V12, V11_inv_V12)

  # Conditional X2: X2_tilde_cond = X2_tilde - V21 V11^{-1} X1_tilde
  V21 <- t(V12)
  X2_tilde_cond <- X2_tilde - V21 %*% V11_inv_X1

  # Stage 2 conditional: I2_eff = X2_tilde_cond' S^{-1} X2_tilde_cond
  R_S <- chol(S)
  S_inv_X2 <- backsolve(R_S, backsolve(R_S, X2_tilde_cond, transpose = TRUE))
  I2_eff <- as.numeric(crossprod(X2_tilde_cond, S_inv_X2))

  I_eff <- I1_eff + I2_eff

  list(
    I1_eff = I1_eff,
    I2_eff = I2_eff,
    I_eff = I_eff,
    X_trt_tilde = X_trt_tilde  # Return projected treatment for verification
  )
}

#' Calibrate the stage 2 boundary given a stage 1 efficacy boundary
#'
#' Solves, under H0,
#'   P(|Z1| > c1) + int_R P(|w1 z1 + w2 Z_{2|1}| > c2) phi(z1) dz1 = alpha
#' for c2, where R is the region over which the trial may continue.
#'
#' @param w1_ref Pre-planned stage 1 combination weight
#' @param efficacy_boundary Stage 1 efficacy boundary c1
#' @param alpha Two-sided type I error rate
#' @param continuation NULL for non-binding futility, in which case R is the
#'   whole of [-c1, c1]. For binding futility, a list of numeric length-2
#'   vectors giving the continuation intervals, as returned by
#'   `continuation_intervals()`.
#' @return Scalar c2
#' @export
calibrate_c2 <- function(w1_ref, efficacy_boundary, alpha = 0.05,
                         continuation = NULL) {

  w2 <- sqrt(1 - w1_ref^2)
  c1 <- efficacy_boundary
  spend1 <- 2 * (1 - pnorm(c1))

  if (spend1 >= alpha) {
    stop(sprintf(
      "Stage 1 boundary c1 = %.3f spends %.4f, exceeding alpha = %.3f. %s",
      c1, spend1, alpha,
      "A larger c1 (smaller w1) or larger alpha is required."))
  }

  cte <- function(z, c2) {
    pnorm((-c2 - w1_ref * z) / w2) + 1 - pnorm((c2 - w1_ref * z) / w2)
  }

  regions <- if (is.null(continuation)) list(c(-c1, c1)) else continuation
  regions <- Filter(function(r) r[2] > r[1], regions)
  if (!length(regions))
    stop("continuation region is empty; alpha cannot be attained at stage 2")

  stage2 <- function(c2) {
    sum(vapply(regions, function(r)
      integrate(function(z) cte(z, c2) * dnorm(z), r[1], r[2])$value,
      numeric(1)))
  }

  ## With binding futility the solution can fall below z_{alpha/2}: the
  ## futility boundary returns more alpha than the interim efficacy stop
  ## consumes. The search interval must allow for that.
  uniroot(function(c2) spend1 + stage2(c2) - alpha, c(0.2, 8))$root
}

#' Continuation intervals from a decision-rule grid
#'
#' Converts a logical `continue` indicator over a z1 grid into a list of
#' intervals, taking interval edges at the midpoints between grid nodes.
#' Returns a list of length-2 numeric vectors; more than one element indicates
#' a non-contiguous continuation region.
#' @export
continuation_intervals <- function(z1_grid, continue_flag) {
  o <- order(z1_grid)
  z <- z1_grid[o]; f <- as.logical(continue_flag)[o]
  if (!any(f)) return(list())

  r <- rle(f)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  lapply(which(r$values), function(i) {
    s <- starts[i]; e <- ends[i]
    lo <- if (s == 1L)        z[1]         else mean(z[c(s - 1L, s)])
    hi <- if (e == length(z)) z[length(z)] else mean(z[c(e, e + 1L)])
    c(lo, hi)
  })
}

# =============================================================================
# Model Interface Specification
# =============================================================================

#' Model functions must return a list with at minimum:
#'   - I1_eff: Stage 1 efficient information (scalar)
#'   - I2_eff: Stage 2 conditional efficient information (scalar)
#'   - I_eff: Total efficient information (scalar)
#'   - w1, w2: Combination weights
#'   - b1: Treatment effect (delta) for power calculations
#'
#' Additional fields can be included for:
#'   - Resource tracking (n1, n2, k1, k2, etc.)
#'   - Cost calculation (any fields used by cost_fn)
#'   - Sample size computation (any fields used by sample_size_fn)
#'
#' The function signature should be: function(design_params, fixed_params)
#' where design_params is a single row from the design grid.

# =============================================================================
# Core Computational Functions
# =============================================================================

#' #' Conditional power (two-sided)
#' #' @param z1 Stage 1 z-statistic
#' #' @param I2_eff Stage 2 conditional efficient information
#' #' @param w1, w2 Combination weights
#' #' @param delta Effect size on linear predictor scale
#' conditional_power <- function(z1, I2_eff, w1, w2, delta) {
#'   mu <- delta * sqrt(I2_eff)
#'   upper_thresh <- (1.96 - w1 * z1) / w2
#'   lower_thresh <- (-1.96 - w1 * z1) / w2
#'   pnorm(mu - upper_thresh) + pnorm(lower_thresh - mu)
#' }

#' Vectorised conditional power over z1 and design grid
#' @param z1_vec Vector of z1 values
#' @param model_summaries Data frame with I2_eff, w1, w2, b1 columns
#' @return Matrix: rows = z1 values, cols = designs
conditional_power_matrix <- function(z1_vec, model_summaries, c2 = NULL,
                                     cp_direction = c("target", "both")) {
  cp_direction <- match.arg(cp_direction)
  nz <- length(z1_vec)
  ng <- nrow(model_summaries)

  W1 <- matrix(model_summaries$w1, nz, ng, byrow = TRUE)
  W2 <- matrix(model_summaries$w2, nz, ng, byrow = TRUE)
  Z1 <- matrix(z1_vec, nz, ng)

  # Non-centrality parameter for z2|1
  NCP2 <- matrix(model_summaries$b1 * sqrt(model_summaries$I2_eff), nz, ng, byrow = TRUE)

  # Degrees of freedom
  df_full <- model_summaries$df_full
  if (is.null(df_full)) df_full <- rep(Inf, ng)

  # df for stage 2 conditional statistic
  # This is the df that would be used for testing with stage 2 data
  df_s2 <- model_summaries$df_s2
  if (is.null(df_s2)) df_s2 <- df_full

  # Critical values from full data df
  if (is.null(c2)) c2 <- qnorm(0.975)
  C2 <- matrix(c2, nz, ng, byrow = TRUE)

  # Thresholds for Z_{2|1} on the z scale
  upper_z <- ( C2 - W1 * Z1) / W2
  lower_z <- (-C2 - W1 * Z1) / W2

  # Section 3.1: Z_{2|1} = Phi^{-1}(F_{t_nu}(T_{2|1})), so the event
  # {Z_{2|1} > u} is the event {T_{2|1} > qt(Phi(u), nu)}. Map the z-scale
  # thresholds onto the t scale before evaluating the non-central t.
  cp_mat <- matrix(NA, nz, ng)

  for (j in seq_len(ng)) {
    ncp_j <- NCP2[1, j]
    df_j  <- df_s2[j]

    if (is.finite(df_j)) {
      upper_t <- qt(pmin(pmax(pnorm(upper_z[, j]), 1e-15), 1 - 1e-15), df = df_j)
      lower_t <- qt(pmin(pmax(pnorm(lower_z[, j]), 1e-15), 1 - 1e-15), df = df_j)
    } else {
      upper_t <- upper_z[, j]
      lower_t <- lower_z[, j]
    }

    if (cp_direction == "target") {
      s <- sign(ncp_j)
      if (s == 0) s <- 1
      upper_s <- if (s > 0) upper_t else -lower_t
      cp_mat[, j] <- suppressWarnings(
        pt(upper_s, df = df_j, ncp = abs(ncp_j), lower.tail = FALSE))
    } else {
      cp_mat[, j] <- suppressWarnings(
        pt(lower_t, df = df_j, ncp = ncp_j) +
          pt(upper_t, df = df_j, ncp = ncp_j, lower.tail = FALSE))
    }
  }

  cp_mat
}

# =============================================================================
# Design Grid and Model Precomputation
# =============================================================================

#' Precompute models for all designs in grid
#' @param design_grid Data frame of design parameters to search over
#' @param model_fn Model function with signature function(design_row, fixed_params)
#' @param fixed_params List of fixed parameters passed to model_fn
#' @param parallel Use parallel computation
#' @param n_cores Number of cores (default: all - 1)
#' @return List of model outputs
precompute_models <- function(design_grid, model_fn, fixed_params,
                              parallel = FALSE, n_cores = NULL, verbose = FALSE) {
  compute_one <- function(i) {
    model_fn(design_grid[i, , drop = FALSE], fixed_params)
  }
  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    model_list <- parallel::mclapply(1:nrow(design_grid), compute_one,
                                     mc.cores = n_cores)
  } else {
    # Use pbapply if available for progress bar
    if (requireNamespace("pbapply", quietly = TRUE) & verbose) {
      model_list <- pbapply::pblapply(1:nrow(design_grid), compute_one)
    } else {
      model_list <- lapply(1:nrow(design_grid), compute_one)
    }
  }

  model_list
}

#' Extract summary statistics from model list
#'
#' @param model_list Output from precompute_models
#' @param cost_fn Cost function: function(model_output, ...) -> scalar RAW cost (unscaled)
#' @param cost_params Fixed cost parameters (e.g., list(rho = 25)), NOT including lambda
#' @param lambda Current lambda value for scaling
#' @param resource_vars Character vector of resource variable names to extract
#' @return Data frame with I2_eff, w1, w2, raw_cost, cost, b1, and any resource_vars
summarise_models <- function(model_list, cost_fn, cost_params, lambda = 1,
                             resource_vars = NULL) {
  # Compute raw (unscaled) costs
  raw_costs <- sapply(model_list, function(d) do.call(cost_fn, c(list(d), cost_params)))

  # Core required fields
  base_df <- data.frame(
    I2_eff = sapply(model_list, function(d) d$I2_eff),
    w1 = sapply(model_list, function(d) d$w1),
    w2 = sapply(model_list, function(d) d$w2),
    raw_cost = raw_costs,
    cost = lambda * raw_costs,
    b1 = sapply(model_list, function(d) d$b1),
    df_s1 = sapply(model_list, function(d) d$df_s1 %||% Inf),
    df_full = sapply(model_list, function(d) d$df_full %||% Inf)
  )

  # Add resource variables if specified
  if (!is.null(resource_vars)) {
    for (var in resource_vars) {
      base_df[[var]] <- sapply(model_list, function(d) {
        if (var %in% names(d)) d[[var]] else NA
      })
    }
  }

  base_df
}

# =============================================================================
# Optimal Design Selection
# =============================================================================

#' Find optimal stage 2 design for each z1 value
#'
#' @param z1_vec Vector of z1 values
#' @param design_grid Design parameter grid
#' @param model_summaries Output from summarise_models
#' @param resource_vars Character vector of resource variable names to include in output
#' @return Data frame with z1, optimal design index, cp, continue flag, and resources
find_optimal_designs <- function(z1_vec, design_grid, model_summaries,
                                 resource_vars = NULL,
                                 w1_ref = NULL,
                                 cost_cap = NULL,
                                 c2 = NULL,
                                 efficacy_boundary = NULL,
                                 cp_direction = NULL) {

  cp_mat <- conditional_power_matrix(z1_vec, model_summaries, c2 = c2,
                                     cp_direction = cp_direction)

  if (is.null(cost_cap)) {
    # Original: λ-penalised criterion
    criterion_mat <- cp_mat - matrix(model_summaries$cost,
                                     nrow = length(z1_vec),
                                     ncol = nrow(design_grid),
                                     byrow = TRUE)
    best_idx <- max.col(criterion_mat)
    best_criterion <- criterion_mat[cbind(seq_along(z1_vec), best_idx)]
    best_cp <- cp_mat[cbind(seq_along(z1_vec), best_idx)]
    futility_criterion <- best_criterion < 0

  } else {
    # Budget-constrained: max CP subject to raw_cost <= cost_cap
    feasible <- model_summaries$raw_cost <= cost_cap

    if (!any(feasible)) {
      # No feasible designs — everything stops
      best_idx <- rep(1L, length(z1_vec))
      best_cp <- rep(0, length(z1_vec))
      best_criterion <- rep(-Inf, length(z1_vec))
      futility_criterion <- rep(TRUE, length(z1_vec))
    } else {
      cp_mat_feas <- cp_mat
      cp_mat_feas[, !feasible] <- -Inf
      best_idx <- max.col(cp_mat_feas)
      best_cp <- cp_mat[cbind(seq_along(z1_vec), best_idx)]
      best_criterion <- best_cp
      # No explicit futility — budget constraint does the work
      # But stop if CP is negligible (no point continuing)
      futility_criterion <- best_cp < 1e-6
    }
  }

  # --- everything below here is unchanged ---
  if (is.null(w1_ref)) {
    w1_ref <- model_summaries$w1[1]
  }

  if (is.null(efficacy_boundary)) {
    efficacy_boundary <- qt(0.975, df = model_summaries$df_full[1]) / w1_ref
  }
  efficacy_stop <- abs(z1_vec) > efficacy_boundary
  futility_stop <- !efficacy_stop & futility_criterion
  continue <- !efficacy_stop & !futility_stop

  result <- data.frame(
    z1 = z1_vec,
    efficacy_stop = efficacy_stop,
    futility_stop = futility_stop,
    continue = continue,
    cp = ifelse(continue, best_cp, ifelse(efficacy_stop, 1, 0)),
    criterion = best_criterion,
    design_idx = ifelse(continue, best_idx, NA),
    w1_ref = w1_ref,
    efficacy_boundary = efficacy_boundary,
    c2 = if (is.null(c2)) NA_real_ else c2
  )

  design_cols <- names(design_grid)
  for (col in design_cols) {
    result[[col]] <- ifelse(continue, design_grid[[col]][best_idx], NA)
  }

  if (!is.null(resource_vars)) {
    for (var in resource_vars) {
      if (var %in% names(model_summaries)) {
        result[[var]] <- ifelse(continue, model_summaries[[var]][best_idx],
                                if (var %in% c("n2", "k2", "m2")) 0 else NA)
      }
    }
  }

  result
}

# =============================================================================
# Power and Expected Sample Size
# =============================================================================

#' Compute power and expected resources using Gauss-Legendre quadrature
#'
#' @param design_grid Design grid
#' @param model_summaries Model summaries
#' @param stage1_info I1_eff from stage 1
#' @param resource_vars Variables to compute expectations for
#' @param z1_range Integration range for z1 (default c(-4, 4))
#' @param n_quad Number of quadrature points
#' @return List with power components, probabilities, and expected resources
compute_power_and_expectations <- function(design_grid, model_summaries, I1_eff,
                                           resource_vars = NULL,
                                           z1_range = NULL,
                                           n_quad = 100,
                                           w1_ref = NULL,
                                           df_s1 = Inf,
                                           df_full = Inf,
                                           cost_cap = NULL,
                                           c2 = NULL,
                                           cp_direction = NULL) {

  b1 <- model_summaries$b1[1]
  mu1 <- b1 * sqrt(I1_eff)

  # Centre the quadrature on mu1. A fixed c(-4, 4) truncates ~5% of the mass
  # when mu1 is around 2.3, which does not affect power_total (power_stage1 is
  # analytic and the continuation region lies inside the range) but does make
  # prob_efficacy + prob_continue + prob_futility sum to less than 1.
  if (is.null(z1_range)) z1_range <- mu1 + c(-8, 8)
  # Gauss-Legendre quadrature
  gl <- statmod::gauss.quad(n_quad, kind = "legendre")
  z1_grid <- (gl$nodes + 1) / 2 * diff(z1_range) + z1_range[1]
  gl_weights <- gl$weights / 2 * diff(z1_range)
  dens    <- z1_density(z1_grid, mu1, df_s1)     # was dnorm(z1_grid - mu1)
  weights <- gl_weights * dens

  # Get w1 if not provided
  if (is.null(w1_ref)) {
    stop("w1_ref must be supplied")
  }

  # Get df if not provided
  if (is.null(df_s1) || is.infinite(df_s1)) {
    df_s1 <- model_summaries$df_s1[1] %||% Inf
  }
  if (is.null(df_full) || is.infinite(df_full)) {
    df_full <- model_summaries$df_full[1] %||% Inf
  }

  t_crit <- qnorm(0.975)
  efficacy_boundary <- t_crit / w1_ref

  # Stage 2 boundary c2, calibrated jointly with c1
  if (is.null(c2)) c2 <- calibrate_c2(w1_ref, efficacy_boundary, alpha = 0.05)
  opt_designs <- find_optimal_designs(z1_grid, design_grid, model_summaries,
                                      resource_vars, w1_ref = w1_ref,
                                      cost_cap = cost_cap,
                                      c2 = c2,
                                      efficacy_boundary = efficacy_boundary,
                                      cp_direction = cp_direction)

  # Stage 1 power: P(|t1| > efficacy_boundary) under H1
  # t1 ~ t(df_s1, ncp = mu1) approximately, or using normal approximation:
  # z1 ~ N(mu1, 1) and we reject if |z1| > efficacy_boundary


  power_stage1 <- stage1_power(efficacy_boundary, mu1, df_s1)

  # Stage 2 power: CP integrated over continuation region
  continue <- opt_designs$continue
  power_stage2 <- sum(opt_designs$cp[continue] * weights[continue])

  power_total <- power_stage1 + power_stage2

  # Probability of each outcome
  prob_efficacy <- sum(weights[opt_designs$efficacy_stop])
  prob_futility <- sum(weights[opt_designs$futility_stop])
  prob_continue <- sum(weights[continue])

  prob_total <- prob_efficacy + prob_futility + prob_continue
  if (abs(prob_total - 1) > 1e-3) {
    warning(sprintf(
      "decision probabilities sum to %.4f; widen z1_range or raise n_quad",
      prob_total))
  }

  # Expected values over continuation region
  expected <- list()
  if (any(continue)) {
    for (col in names(opt_designs)) {
      if (is.numeric(opt_designs[[col]]) && !col %in% c("z1", "efficacy_stop",
                                                        "futility_stop", "continue")) {
        vals <- opt_designs[[col]]
        expected[[paste0("E_", col)]] <- sum(vals[continue] * weights[continue], na.rm = TRUE) /
          sum(weights[continue])
      }
    }
  }

  list(
    power = list(
      total = power_total,
      stage1 = power_stage1,
      stage2 = power_stage2
    ),
    probabilities = list(
      efficacy_stop = prob_efficacy,
      futility_stop = prob_futility,
      continue = prob_continue
    ),
    expected = expected,
    quadrature = list(
      z1 = z1_grid,
      weights = weights,
      gl_weights_raw = gl_weights,
      optimal_designs = opt_designs
    ),
    params = list(
      mu1 = mu1,
      w1_ref = w1_ref,
      efficacy_boundary = efficacy_boundary,
      t_crit = t_crit,
      c2 = c2,
      df_s1 = df_s1,
      df_full = df_full
    )
  )
}

# =============================================================================
# Lambda Calibration (Bisection Search)
# =============================================================================

#' Find lambda that achieves target power
#'
#' @param target_power Target overall power (e.g., 0.8)
#' @param design_grid Data frame of stage 2 designs
#' @param model_fn Model function
#' @param fixed_params Fixed parameters for model
#' @param cost_fn Cost function: function(model_output, ...) -> raw cost (unscaled)
#' @param cost_params Fixed cost parameters (e.g., list(rho = 25)), NOT including lambda
#' @param resource_vars Variables to track expectations for
#' @param lambda_range Initial search range
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
#' @param z1_range Range for z1 integration
#' @param n_quad Quadrature points
#' @param verbose Print progress
#' @return List with lambda, power, and full results
.find_lambda_inner <- function(target_power = 0.8,
                                  design_grid,
                                  model_fn,
                                  fixed_params,
                                  cost_fn,
                                  cost_params = list(rho = 25),
                                  resource_vars = NULL,
                                  method = c("lambda", "cost_cap"),
                                  lambda_range = c(1e-8, 1),
                                  tol = 0.01,
                                  max_iter = 50,
                                  z1_range = NULL,
                                  n_quad = 100,
                                  verbose = TRUE,
                                  # --- Interim overrides (lock to planned values) ---
                                  w1_override = NULL,
                                  efficacy_override = NULL,
                                  t_crit_override = NULL,
                                  c2_override = NULL,
                                  futility = c("nonbinding", "binding"),
                                  c2_tol = 1e-4,
                                  c2_max_iter = 5,
                                  c2_fixed = NULL) {

  method <- match.arg(method)

  # Precompute models once (expensive part)
  if (verbose) cat("Precomputing models...\n")
  model_list <- precompute_models(design_grid, model_fn, fixed_params,
                                  parallel = FALSE, n_cores = NULL, verbose = verbose)

  I1_all <- vapply(model_list, function(m) m$I1_eff, numeric(1))
  b1_all <- vapply(model_list, function(m) m$b1,     numeric(1))
  if (diff(range(I1_all)) > 1e-2 * mean(I1_all))
    stop(sprintf("I1_eff varies across the stage 2 grid (%.6g to %.6g)",
                 min(I1_all), max(I1_all)))
  if (diff(range(b1_all)) > 1e-10) stop("b1 varies across the stage 2 grid")
  I1_eff <- I1_all[1]; b1 <- b1_all[1]; mu1 <- b1 * sqrt(I1_eff)
  I1_eff <- model_list[[1]]$I1_eff
  I_eff <- model_list[[1]]$I_eff
  b1 <- model_list[[1]]$b1

  if (!is.null(w1_override)) {
    # Interim mode: use pre-planned weights, recompute mu1 from updated I1
    w1_ref <- w1_override
    t_crit_s1 <- t_crit_override %||% 1.96
    efficacy_boundary <- efficacy_override %||% (t_crit_s1 / w1_ref)
    mu1 <- b1 * sqrt(I1_eff)  # updated I1_eff with new ICC
    df_s1 <- Inf  # not used when overriding
  } else {
    df_s1     <- model_list[[1]]$df_s1 %||% Inf
    t_crit_s1 <- qnorm(0.975)                        # see (3)
    I2_all    <- vapply(model_list, function(m) m$I2_eff, numeric(1))
    ms0       <- summarise_models(model_list, cost_fn, cost_params,
                                  lambda_range[1], resource_vars)
    eval_w1 <- function(w) compute_power_and_expectations(
      design_grid, ms0, I1_eff, resource_vars, z1_range, n_quad,
      w1_ref = w, df_s1 = df_s1, df_full = NULL,
      c2 = calibrate_c2(w, t_crit_s1 / w, alpha = 0.05))$power$total
    lo <- sqrt(I1_eff / (I1_eff + max(I2_all)))
    hi <- sqrt(I1_eff / (I1_eff + min(I2_all)))
    w1_ref <- if (isTRUE(all.equal(lo, hi))) lo
    else optimize(function(w) -eval_w1(w), c(lo, hi))$minimum
    efficacy_boundary <- t_crit_s1 / w1_ref
  }
  futility <- match.arg(futility)
  c2 <- c2_fixed %||% c2_override %||%
    calibrate_c2(w1_ref, efficacy_boundary, alpha = 0.05)

  t_crit_s1 <- qnorm(0.975)                       # was qt(0.975, df = df_s1)
  efficacy_boundary <- t_crit_s1 / w1_ref
  efficacy_boundary_t <- if (is.finite(df_s1))    # what the trial's t stat is compared to
    qt(pnorm(efficacy_boundary), df = df_s1) else efficacy_boundary
  power_stage1 <- stage1_power(efficacy_boundary, mu1, df_s1)

  if (verbose) {
    cat(sprintf("Method: %s\n", method))
    cat(sprintf("df_s1 = %.0f, t_crit = %.3f\n", df_s1, t_crit_s1))
    cat(sprintf("w1 = %.4f, efficacy boundary |z1| > %.3f\n", w1_ref, efficacy_boundary))
    cat(sprintf("stage 2 boundary |Z| > %.3f (calibrated; z_a/2 = %.3f)\n",
                c2, qnorm(0.975)))
    cat(sprintf("Stage 1 power: %.4f\n", power_stage1))
  }

  # Precompute model summaries (lambda=1 as placeholder; raw_cost is what matters)
  model_summaries <- summarise_models(model_list, cost_fn, cost_params,
                                      lambda = 1, resource_vars)

  build_results <- function(lambda, model_summaries, power_results) {
    od <- power_results$quadrature$optimal_designs
    zg <- power_results$quadrature$z1
    stopifnot(!is.null(od$futility_stop), length(zg) == length(od$futility_stop))
    fz <- zg[od$futility_stop & zg < efficacy_boundary]
    futility_boundary <- if (length(fz)) max(fz) else -Inf
    if (!is.finite(futility_boundary))
      warning("no futility region found on the quadrature grid; ",
              "protocol futility rule will never fire")
    list(
      power = power_results$power,
      probabilities = power_results$probabilities,
      expected = power_results$expected,
      quadrature = power_results$quadrature,
      models = list(
        list = model_list,
        summaries = model_summaries,
        design_grid = design_grid
      ),
      params = list(
        mu1 = mu1,
        b1 = b1,
        I1_eff = I1_eff,
        I_eff = I_eff,
        w1_ref = w1_ref,
        df_s1 = df_s1,
        t_crit_s1 = t_crit_s1,
        efficacy_boundary = efficacy_boundary,
        efficacy_boundary_t = efficacy_boundary_t,
        futility_boundary = futility_boundary,
        c2 = c2,
        futility = futility,
        lambda = lambda,
        cost_params = cost_params,
        method = method
      )
    )
  }

  if (power_stage1 >= target_power) {
    if (verbose) cat("Target power already achieved by stage 1 alone.\n")
    return(list(
      converged = TRUE,
      lambda = Inf,
      cost_cap = 0,
      cost_params = cost_params,
      power = power_stage1,
      target = target_power,
      method = method,
      futility = futility,
      c2 = c2,
      message = "Target achieved by stage 1 efficacy stopping alone",
      results = NULL,
      iterations = 0
    ))
  }

  # --- Dispatch by method ---

  if (method == "lambda") {
    # ---- Existing λ-bisection (unchanged logic) ----

    compute_for_lambda <- function(lambda) {
      ms <- summarise_models(model_list, cost_fn, cost_params, lambda, resource_vars)
      pr <- compute_power_and_expectations(
        design_grid, ms, I1_eff, resource_vars,
        z1_range, n_quad,
        w1_ref = w1_ref, df_s1 = df_s1,
        df_full = ms$df_full[1],
        c2 = c2
      )
      list(power = pr$power$total, summaries = ms, results = pr)
    }

    if (verbose) cat("Checking maximum achievable power (lambda -> 0)...\n")
    res_max <- compute_for_lambda(lambda_range[1])
    power_max <- res_max$power
    if (verbose) cat(sprintf("Maximum achievable power: %.4f\n", power_max))

    if (power_max < target_power) {
      if (verbose) cat(sprintf("Target power %.3f not achievable. Maximum is %.4f\n",
                               target_power, power_max))
      return(list(
        converged = FALSE, lambda = lambda_range[1], cost_cap = NA,
        cost_params = cost_params, power = power_max, target = target_power,
        method = method,
        message = sprintf("Target power not achievable. Maximum power is %.4f", power_max),
        results = build_results(lambda_range[1], res_max$summaries, res_max$results),
        iterations = 1
      ))
    }

    res_min <- compute_for_lambda(lambda_range[2])
    power_min <- res_min$power

    if (power_min > target_power) {
      if (verbose) cat("Expanding search range upward...\n")
      while (power_min > target_power && lambda_range[2] < 1e6) {
        lambda_range[2] <- lambda_range[2] * 10
        res_min <- compute_for_lambda(lambda_range[2])
        power_min <- res_min$power
      }
    }

    if (verbose) cat(sprintf("Power range: [%.4f, %.4f] for lambda in [%.2e, %.2e]\n",
                             power_min, power_max, lambda_range[1], lambda_range[2]))

    lambda_lo <- lambda_range[1]
    lambda_hi <- lambda_range[2]
    history <- data.frame(iteration = integer(), lambda = numeric(), power = numeric())
    res_mid <- NULL; lambda_mid <- NULL; power_mid <- NULL

    for (iter in 1:max_iter) {
      lambda_mid <- exp((log(lambda_lo) + log(lambda_hi)) / 2)
      res_mid <- compute_for_lambda(lambda_mid)
      power_mid <- res_mid$power

      history <- rbind(history, data.frame(iteration = iter, lambda = lambda_mid,
                                           power = power_mid))
      if (verbose) cat(sprintf("Iter %2d: lambda = %.4e, power = %.4f (target = %.4f)\n",
                               iter, lambda_mid, power_mid, target_power))

      if (abs(power_mid - target_power) < tol) {
        if (verbose) cat("Converged!\n")
        return(list(
          converged = TRUE, lambda = lambda_mid, cost_cap = NA,
          cost_params = cost_params, power = power_mid, target = target_power,
          method = method, message = "Converged",
          results = build_results(lambda_mid, res_mid$summaries, res_mid$results),
          iterations = iter, history = history
        ))
      }

      if (power_mid > target_power) lambda_lo <- lambda_mid else lambda_hi <- lambda_mid
    }

    warning("Maximum iterations reached without convergence")
    list(
      converged = FALSE, lambda = lambda_mid, cost_cap = NA,
      cost_params = cost_params, power = power_mid, target = target_power,
      method = method,
      message = sprintf("Max iterations reached. Power = %.4f", power_mid),
      results = build_results(lambda_mid, res_mid$summaries, res_mid$results),
      iterations = max_iter, history = history
    )

  } else {
    # ---- Budget-constrained bisection ----

    raw_costs <- model_summaries$raw_cost
    cost_range <- range(raw_costs)

    compute_for_cap <- function(cap) {
      pr <- compute_power_and_expectations(
        design_grid, model_summaries, I1_eff, resource_vars,
        z1_range, n_quad,
        w1_ref = w1_ref, df_s1 = df_s1,
        df_full = model_summaries$df_full[1],
        cost_cap = cap,
        c2 = c2
      )
      list(power = pr$power$total, summaries = model_summaries, results = pr)
    }

    # Check max power (full budget)
    if (verbose) cat("Checking maximum achievable power (full budget)...\n")
    res_max <- compute_for_cap(cost_range[2])
    power_max <- res_max$power
    if (verbose) cat(sprintf("Maximum achievable power: %.4f\n", power_max))

    if (power_max < target_power) {
      if (verbose) cat(sprintf("Target power %.3f not achievable. Maximum is %.4f\n",
                               target_power, power_max))
      return(list(
        converged = FALSE, lambda = NA, cost_cap = cost_range[2],
        cost_params = cost_params, power = power_max, target = target_power,
        method = method,
        message = sprintf("Target power not achievable. Maximum power is %.4f", power_max),
        results = build_results(NA, model_summaries, res_max$results),
        iterations = 1
      ))
    }

    # Check min power (cheapest design only)
    res_min <- compute_for_cap(cost_range[1])
    power_min <- res_min$power

    if (verbose) cat(sprintf("Power range: [%.4f, %.4f] for cost_cap in [%.1f, %.1f]\n",
                             power_min, power_max, cost_range[1], cost_range[2]))

    cap_lo <- cost_range[1]
    cap_hi <- cost_range[2]
    history <- data.frame(iteration = integer(), cost_cap = numeric(), power = numeric())
    res_mid <- NULL; cap_mid <- NULL; power_mid <- NULL

    for (iter in 1:max_iter) {
      cap_mid <- (cap_lo + cap_hi) / 2
      res_mid <- compute_for_cap(cap_mid)
      power_mid <- res_mid$power

      history <- rbind(history, data.frame(iteration = iter, cost_cap = cap_mid,
                                           power = power_mid))
      if (verbose) cat(sprintf("Iter %2d: cost_cap = %.1f, power = %.4f (target = %.4f)\n",
                               iter, cap_mid, power_mid, target_power))

      if (abs(power_mid - target_power) < tol) {
        if (verbose) cat("Converged!\n")
        return(list(
          converged = TRUE, lambda = NA, cost_cap = cap_mid,
          cost_params = cost_params, power = power_mid, target = target_power,
          method = method, message = "Converged",
          results = build_results(NA, model_summaries, res_mid$results),
          iterations = iter, history = history
        ))
      }

      if (power_mid > target_power) cap_hi <- cap_mid else cap_lo <- cap_mid
    }

    warning("Maximum iterations reached without convergence")
    list(
      converged = FALSE, lambda = NA, cost_cap = cap_mid,
      cost_params = cost_params, power = power_mid, target = target_power,
      method = method,
      message = sprintf("Max iterations reached. Power = %.4f", power_mid),
      results = build_results(NA, model_summaries, res_mid$results),
      iterations = max_iter, history = history
    )
  }
}

#' Find lambda that achieves target power
#'
#' Wrapper around the solve that resolves the circular dependence between the
#' stage 2 boundary and the futility region when binding futility is used:
#' c2 depends on the continuation region, which depends on lambda, which is
#' calibrated against a power target computed using c2. Converges quickly
#' because the futility boundary moves only weakly with c2.
#'
#' @inheritParams .find_lambda_inner
#' @export
find_lambda_for_power <- function(..., futility = c("nonbinding", "binding"),
                                  c2_tol = 1e-4, c2_max_iter = 5,
                                  verbose = FALSE) {

  futility <- match.arg(futility)

  ## Non-binding: single pass, identical to previous behaviour.
  if (futility == "nonbinding") {
    return(.find_lambda_inner(..., futility = "nonbinding", verbose = verbose))
  }

  fit <- NULL
  c2  <- NULL

  for (it in seq_len(c2_max_iter)) {

    fit <- .find_lambda_inner(..., futility = "binding", verbose = verbose,
                              c2_fixed = c2)

    ## No continuation region to read (e.g. stage 1 alone reaches the target)
    if (is.null(fit$results)) return(fit)

    p <- fit$results$params

    ## Calibrate against the region interim_analysis() actually applies: the
    ## single interval (c_f, c1). Using the raw `continue` flags instead can
    ## span a wider, non-contiguous region on the quadrature grid, giving a c2
    ## that is too low and an attained error rate above alpha.
    if (!is.finite(p$futility_boundary) ||
        p$futility_boundary >= p$efficacy_boundary) {
      warning("no valid futility boundary; retaining non-binding c2")
      return(fit)
    }
    cont <- list(c(p$futility_boundary, p$efficacy_boundary))

    c2_new <- calibrate_c2(p$w1_ref, p$efficacy_boundary, alpha = 0.05,
                           continuation = cont)

    if (verbose) {
      cat(sprintf("c2 iteration %d: %s -> %.4f  continuation %s\n", it,
                  if (is.null(c2)) "(non-binding start)" else sprintf("%.4f", c2),
                  c2_new,
                  paste(sprintf("[%.2f, %.2f]",
                                vapply(cont, `[`, numeric(1), 1),
                                vapply(cont, `[`, numeric(1), 2)),
                        collapse = " ")))
    }

    if (!is.null(c2) && abs(c2_new - c2) < c2_tol) {
      c2 <- c2_new
      break
    }
    c2 <- c2_new
  }

  ## Final pass so the returned object is solved at the converged c2
  fit <- .find_lambda_inner(..., futility = "binding", verbose = verbose,
                            c2_fixed = c2)
  if (!isTRUE(all.equal(fit$results$params$c2, c2)))
    stop(sprintf("binding c2 not propagated: params$c2 = %s, expected %.4f",
                 format(fit$results$params$c2), c2))
  fit$c2_iterations <- it
  fit
}

# =============================================================================
# Sample Size Computation (Generic)
# =============================================================================

#' Compute sample size distribution and summaries
#'
#' @param results Output from find_lambda_for_power or adaptive_trial_analysis
#' @param sample_size_fn Function that computes sample sizes.
#'        Signature: function(opt_designs, stage1_info, ...) -> list with:
#'          - n_stage1: Scalar or vector of stage 1 sample sizes
#'          - n_stage2: Vector of stage 2 sample sizes for each z1 point
#'          - metrics: Named list of additional metrics to report (e.g., max_n2, E_clusters)
#' @param hypothesis "H1" (default) or "H0" - determines z1 distribution
#' @param ... Additional arguments passed to sample_size_fn
compute_sample_size <- function(results, sample_size_fn = NULL,
                                hypothesis = c("H1", "H0"), ...) {

  hypothesis <- match.arg(hypothesis)

  if (is.null(sample_size_fn)) {
    stop("sample_size_fn must be provided")
  }

  # Check structure
  if (is.null(results$quadrature)) {
    stop("results$quadrature is NULL - check results structure")
  }
  if (is.null(results$quadrature$optimal_designs)) {
    stop("results$quadrature$optimal_designs is NULL - check compute_power_and_expectations")
  }

  opt <- results$quadrature$optimal_designs
  z1_grid <- results$quadrature$z1

  if (nrow(opt) == 0) {
    stop("optimal_designs has 0 rows")
  }
  if (length(z1_grid) != nrow(opt)) {
    stop(sprintf("Mismatch: z1_grid has %d elements but optimal_designs has %d rows",
                 length(z1_grid), nrow(opt)))
  }

  # Recompute weights under specified hypothesis
  mu1 <- if (hypothesis == "H1") results$params$mu1 else 0

  # Get raw quadrature weights
  gl_weights_raw <- results$quadrature$gl_weights_raw
  if (is.null(gl_weights_raw)) {
    # Fallback: estimate from stored weights
    gl_weights_raw <- results$quadrature$weights / dnorm(z1_grid - results$params$mu1)
  }
  weights <- gl_weights_raw * dnorm(z1_grid - mu1)
  # In compute_sample_size, before calling sample_size_fn:
  cost_params <- results$params$cost_params %||% list(rho = 30)
  ss_info <- do.call(sample_size_fn, c(list(opt, results), cost_params))

  n_stage1 <- ss_info$n_stage1
  n_stage2 <- ss_info$n_stage2
  n_total <- n_stage1 + n_stage2

  # Ensure vectors
  if (length(n_stage1) == 1) n_stage1 <- rep(n_stage1, length(z1_grid))

  # Distribution data frame
  distribution <- data.frame(
    z1 = z1_grid,
    weight = weights,
    efficacy_stop = opt$efficacy_stop,
    futility_stop = opt$futility_stop,
    continue = opt$continue,
    n_stage1 = n_stage1,
    n_stage2 = n_stage2,
    n_total = n_total
  )

  # Add design columns from opt
  design_cols <- setdiff(names(opt), c("z1", "efficacy_stop", "futility_stop",
                                       "continue", "cp", "criterion", "design_idx"))
  for (col in design_cols) {
    if (col %in% names(opt)) {
      distribution[[col]] <- opt[[col]]
    }
  }

  # Summaries
  E_n_stage2 <- sum(n_stage2 * weights)
  E_n_total <- sum(n_total * weights)

  # Maximum and minimum sample size
  n_total_max <- if (!is.null(ss_info$n_total_max)) ss_info$n_total_max else max(n_total)
  n_total_min <- if (!is.null(ss_info$n_total_min)) ss_info$n_total_min else min(n_total)

  # Probabilities of each outcome
  prob_efficacy <- sum(weights[opt$efficacy_stop])
  prob_futility <- sum(weights[opt$futility_stop])
  prob_continue <- sum(weights[opt$continue])

  # Expected sample size conditional on outcome
  E_n_given_efficacy <- n_stage1[1]
  E_n_given_futility <- n_stage1[1]

  continue_idx <- opt$continue
  E_n_given_continue <- if (any(continue_idx)) {
    sum(n_total[continue_idx] * weights[continue_idx]) / sum(weights[continue_idx])
  } else {
    NA
  }

  # Build summary data frame
  summary_metrics <- c(
    "E_N" = E_n_total,
    "E_N_efficacy" = E_n_given_efficacy,
    "E_N_futility" = E_n_given_futility,
    "E_N_continue" = E_n_given_continue,
    "N_min" = n_total_min,
    "N_max" = n_total_max,
    "N_stage1" = n_stage1[1],
    "E_N_stage2" = E_n_stage2,
    "P_efficacy" = prob_efficacy,
    "P_futility" = prob_futility,
    "P_continue" = prob_continue
  )

  # Add custom metrics from sample_size_fn
  if (!is.null(ss_info$metrics)) {
    summary_metrics <- c(summary_metrics, unlist(ss_info$metrics))
  }

  summary_df <- data.frame(
    metric = names(summary_metrics),
    value = as.numeric(summary_metrics)
  )

  # Quantiles of sample size distribution
  ord <- order(n_total)
  cum_prob <- cumsum(weights[ord]) / sum(weights)
  n_sorted <- n_total[ord]

  quantile_fn <- function(p) {
    idx <- which(cum_prob >= p)[1]
    if (is.na(idx)) n_sorted[length(n_sorted)] else n_sorted[idx]
  }

  list(
    hypothesis = hypothesis,
    mu1 = mu1,
    summary = summary_df,
    quantiles = data.frame(
      p = c(0.1, 0.25, 0.5, 0.75, 0.9),
      n = sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), quantile_fn)
    ),
    distribution = distribution
  )
}

# =============================================================================
# Utility Functions
# =============================================================================

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
