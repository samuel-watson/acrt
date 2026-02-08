# =============================================================================
# Adaptive Cluster Trial Framework v2
# Generalised to support arbitrary stage 2 design parameters
# =============================================================================

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
conditional_power_matrix <- function(z1_vec, model_summaries) {
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
  t_crit <- qt(0.975, df = df_full)
  T_CRIT <- matrix(t_crit, nz, ng, byrow = TRUE)

  # Thresholds for t2|1 to reject
  upper <- (T_CRIT - W1 * Z1) / W2
  lower <- (-T_CRIT - W1 * Z1) / W2

  # Conditional power using non-central t for t2|1
  # t2|1 ~ t(df_s2, ncp = mu2|1)
  cp_mat <- matrix(NA, nz, ng)

  for (j in seq_len(ng)) {
    ncp_j <- NCP2[1, j]
    df_j <- df_s2[j]

    # P(T < lower) + P(T > upper) where T ~ t(df, ncp)
    cp_mat[, j] <- pt(lower[, j], df = df_j, ncp = ncp_j) +
      pt(upper[, j], df = df_j, ncp = ncp_j, lower.tail = FALSE)
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
                              parallel = FALSE, n_cores = NULL) {

  compute_one <- function(i) {
    model_fn(design_grid[i, , drop = FALSE], fixed_params)
  }

  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    model_list <- parallel::mclapply(1:nrow(design_grid), compute_one,
                                     mc.cores = n_cores)
  } else {
    # Use pbapply if available for progress bar
    if (requireNamespace("pbapply", quietly = TRUE)) {
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
                                 cost_cap = NULL) {

  cp_mat <- conditional_power_matrix(z1_vec, model_summaries)

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

  efficacy_boundary <- 1.96 / w1_ref
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
    efficacy_boundary = efficacy_boundary
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
                                           z1_range = c(-4, 4), n_quad = 50,
                                           w1_ref = NULL,
                                           df_s1 = Inf,
                                           df_full = Inf,
                                           cost_cap = NULL) {

  b1 <- model_summaries$b1[1]
  mu1 <- b1 * sqrt(I1_eff)

  # Gauss-Legendre quadrature
  gl <- statmod::gauss.quad(n_quad, kind = "legendre")
  z1_grid <- (gl$nodes + 1) / 2 * diff(z1_range) + z1_range[1]
  gl_weights <- gl$weights / 2 * diff(z1_range)
  dens <- dnorm(z1_grid - mu1)
  weights <- gl_weights * dens

  # Get w1 if not provided
  if (is.null(w1_ref)) {
    w1_ref <- model_summaries$w1[1]
  }

  # Get df if not provided
  if (is.null(df_s1) || is.infinite(df_s1)) {
    df_s1 <- model_summaries$df_s1[1] %||% Inf
  }
  if (is.null(df_full) || is.infinite(df_full)) {
    df_full <- model_summaries$df_full[1] %||% Inf
  }

  # Critical value for efficacy boundary (from full data df)
  t_crit <- qt(0.975, df = df_full)
  efficacy_boundary <- t_crit / w1_ref

  # Find optimal designs with correct boundary
  opt_designs <- find_optimal_designs(z1_grid, design_grid, model_summaries,
                                      resource_vars, w1_ref = w1_ref,
                                      cost_cap = cost_cap)

  # Stage 1 power: P(|t1| > efficacy_boundary) under H1
  # t1 ~ t(df_s1, ncp = mu1) approximately, or using normal approximation:
  # z1 ~ N(mu1, 1) and we reject if |z1| > efficacy_boundary

  if (is.finite(df_s1)) {
    # Use non-central t
    power_stage1 <- pt(-efficacy_boundary, df = df_s1, ncp = mu1) +
      pt(efficacy_boundary, df = df_s1, ncp = mu1, lower.tail = FALSE)
  } else {
    # Fall back to normal
    power_stage1 <- pnorm(-efficacy_boundary - mu1) + (1 - pnorm(efficacy_boundary - mu1))
  }

  # Stage 2 power: CP integrated over continuation region
  continue <- opt_designs$continue
  power_stage2 <- sum(opt_designs$cp[continue] * weights[continue])

  power_total <- power_stage1 + power_stage2

  # Probability of each outcome
  prob_efficacy <- sum(weights[opt_designs$efficacy_stop])
  prob_futility <- sum(weights[opt_designs$futility_stop])
  prob_continue <- sum(weights[continue])

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
find_lambda_for_power <- function(target_power = 0.8,
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
                                  z1_range = c(-4, 4),
                                  n_quad = 50,
                                  verbose = TRUE,
                                  # --- Interim overrides (lock to planned values) ---
                                  w1_override = NULL,
                                  efficacy_override = NULL,
                                  t_crit_override = NULL) {

  method <- match.arg(method)

  # Precompute models once (expensive part)
  if (verbose) cat("Precomputing models...\n")
  model_list <- precompute_models(design_grid, model_fn, fixed_params,
                                  parallel = FALSE, n_cores = NULL)

  I1_eff <- model_list[[1]]$I1_eff
  I_eff <- model_list[[1]]$I_eff
  b1 <- model_list[[1]]$b1
  mu1 <- b1 * sqrt(I1_eff)
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
    # Planning mode: compute everything from scratch
    mu1 <- b1 * sqrt(I1_eff)
    w1_ref <- sqrt(I1_eff / I_eff)
    df_s1 <- model_list[[1]]$df_s1 %||% Inf
    t_crit_s1 <- qt(0.975, df = df_s1)
    efficacy_boundary <- t_crit_s1 / w1_ref
  }

  power_stage1 <- pnorm(-efficacy_boundary - mu1) + (1 - pnorm(efficacy_boundary - mu1))
  if (verbose) {
    cat(sprintf("Method: %s\n", method))
    cat(sprintf("df_s1 = %.0f, t_crit = %.3f\n", df_s1, t_crit_s1))
    cat(sprintf("w1 = %.4f, efficacy boundary |z1| > %.3f\n", w1_ref, efficacy_boundary))
    cat(sprintf("Stage 1 power: %.4f\n", power_stage1))
  }

  # Precompute model summaries (lambda=1 as placeholder; raw_cost is what matters)
  model_summaries <- summarise_models(model_list, cost_fn, cost_params,
                                      lambda = 1, resource_vars)

  build_results <- function(lambda, model_summaries, power_results) {
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
        df_full = ms$df_full[1]
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
        cost_cap = cap
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

  # Call user-provided sample size function
  ss_info <- sample_size_fn(opt, results, ...)

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
