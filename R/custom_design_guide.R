# Custom Design Guide for adaptiveCRT
# ====================================
#
# This guide explains how to create custom adaptive trial designs using the
# adaptiveCRT framework. You only need to write custom code if the built-in
# designs (parallel_crt, crossover_crt, stepped_wedge_crt) don't fit your needs.
#
# For quick reference, run: custom_design_help()
# For a specific topic, run: custom_design_help("model_builder")

# =============================================================================
# OVERVIEW
# =============================================================================
#
# To create a custom design, you need to provide:
#
#   1. model_builder    - Computes efficient information decomposition (REQUIRED)
#   2. Resources list   - Formulas for sample sizes, clusters, etc. (REQUIRED)
#   3. cost_fn          - Only if auto-generated cost doesn't fit your needs
#   4. sample_size_fn   - Only if auto-generated doesn't fit your needs
#
# The framework auto-generates cost and sample size functions from your
# resource formulas, so most users only need to write the model_builder.
#
# =============================================================================


# =============================================================================
# QUICK START TEMPLATE
# =============================================================================

#' Create a custom CRT design (copy and modify this template)
#'
#' @example
#' # 1. Define your model builder (see detailed docs below)
#' my_model_builder <- function(design_params, fixed_params) {
#'   # ... your code here ...
#'   list(I1_eff=..., I2_eff=..., I_eff=..., w1=..., w2=..., b1=...)
#' }
#'
#' # 2. Create design specification
#' my_spec <- crt_design_spec(
#'   stage1_params = c("k1", "m1"),
#'   stage2_params = c("k2", "m2"),
#'   resources = list(
#'     n_s1 = ~ 2 * k1 * m1,
#'     n_s2 = ~ 2 * k1 * m2 + 2 * k2 * m2,
#'     clusters_s1 = ~ 2 * k1,
#'     clusters_s2 = ~ 2 * k2
#'   ),
#'   model_builder = my_model_builder,
#'   design_type = "my_custom_type"
#' )
#'
#' # 3. Create design object
#' my_design <- structure(
#'   list(
#'     spec = my_spec,
#'     fixed_params = list(icc = 0.05, cac = 0.8, delta = 0.25),
#'     stage1_grid = expand.grid(k1 = 8:12, m1 = c(20, 30, 40)),
#'     stage2_grid_fn = function(s1) expand.grid(m2 = 10:50, k2 = 0:4),
#'     rho = 30
#'   ),
#'   class = c("my_custom_crt", "crt_design")
#' )
#'
#' # 4. Analyse
#' results <- adaptive_analysis(my_design, target_power = 0.8)


# =============================================================================
# MODEL_BUILDER FUNCTION
# =============================================================================
#
# The model_builder is the core function that computes the statistical
# quantities needed for adaptive decision-making.
#
# -----------------------------------------------------------------------------
# SIGNATURE
# -----------------------------------------------------------------------------
#
#   model_builder <- function(design_params, fixed_params) { ... }
#
# -----------------------------------------------------------------------------
# INPUTS
# -----------------------------------------------------------------------------
#
# design_params: Single-row data.frame with stage 2 parameters
#   - Comes from your stage2_grid (one row at a time)
#   - Access as: design_params$m2, design_params$k2, etc.
#
# fixed_params: List containing:
#   - Stage 1 parameters (merged from stage1_grid automatically)
#   - Design constants (icc, cac, delta, etc.)
#   - cache: Environment for model caching (IMPORTANT for performance)
#
#   Access as: fixed_params$k1, fixed_params$icc, fixed_params$cache
#
# -----------------------------------------------------------------------------
# REQUIRED OUTPUTS
# -----------------------------------------------------------------------------
#
# Must return a list with these fields:
#
#   I1_eff  : Stage 1 efficient information (scalar > 0)
#   I2_eff  : Stage 2 efficient information (scalar > 0)
#   I_eff   : Total efficient information (scalar > 0)
#   w1      : Weight for stage 1 = sqrt(I1_eff / I_eff)
#   w2      : Weight for stage 2 = sqrt(I2_eff / I_eff)
#   b1      : Treatment effect (same as delta)
#
# Mathematical interpretation:
#   - Z1 ~ N(b1 * sqrt(I1_eff), 1) under H1
#   - Combination test: Z = w1*Z1 + w2*Z2|1
#   - w1^2 + w2^2 = 1 (weights sum to 1 in variance)
#
# -----------------------------------------------------------------------------
# OPTIONAL OUTPUTS (for resource tracking)
# -----------------------------------------------------------------------------
#
# Include any quantities needed by cost or sample size calculations:
#
#   n1, n2      : Observations in each stage
#   m1, m2      : Per-cluster sample sizes
#   k1, k2      : Clusters per arm
#   t1, t2      : Time periods
#   T_total     : Total duration
#   ... any other design-specific values
#
# These are used by:
#   - The auto-generated cost function
#   - The auto-generated sample size function
#   - Your custom functions if you write them
#
# -----------------------------------------------------------------------------
# COMPUTING EFFICIENT INFORMATION
# -----------------------------------------------------------------------------
#
# Use the efficient_score_decomposition() helper:
#
#   # V = marginal covariance matrix of the full design
#   # X = fixed effects design matrix
#   # idx1 = row indices for stage 1
#   # idx2 = row indices for stage 2
#   # j = column index of treatment effect in X
#
#   eff <- efficient_score_decomposition(X, V, idx1, idx2, j = j_treatment)
#
#   # Returns: eff$I1_eff, eff$I2_eff, eff$I_eff
#
# "Efficient" information projects out nuisance parameters (intercepts,
# period effects) that aren't estimable within each stage alone.
#
# -----------------------------------------------------------------------------
# MODEL CACHING (CRITICAL FOR PERFORMANCE)
# -----------------------------------------------------------------------------
#
# The model_builder is called many times (once per stage 2 design option).
# Building glmmrBase models is slow. Use caching:
#
#   cache <- fixed_params$cache
#   cache_key <- paste(k1, k2, sep = "_")  # Unique ID for model structure
#
#   if (is.null(cache$mod) || cache$cache_key != cache_key) {
#     # Model structure changed - must rebuild
#     cache$mod <- Model$new(...)
#     cache$cache_key <- cache_key
#   } else {
#     # Only weights/parameters changed - just update
#     glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, new_weights)
#     cache$mod$update_parameters(cov.pars = new_cov_params)
#   }
#
# Rule of thumb: Rebuild when cluster structure changes, update when only
# sample sizes or variance parameters change.
#
# -----------------------------------------------------------------------------
# COMPUTING V MATRIX FROM glmmrBase MODEL
# -----------------------------------------------------------------------------
#
# Standard pattern:
#
#   X <- cache$mod$mean$X           # Fixed effects design matrix
#   D <- cache$mod$covariance$D     # Random effects covariance
#   Z <- cache$mod$covariance$Z     # Random effects design matrix
#
#   library(Matrix)
#   Z_sp <- Matrix(Z, sparse = TRUE)
#   D_sp <- Matrix(D, sparse = TRUE)
#
#   # Marginal covariance: V = sigma^2 * I + Z D Z'
#   S <- Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
#   V <- as.matrix(S)
#

# =============================================================================
# RESOURCES SPECIFICATION
# =============================================================================
#
# Resources are formulas that define how to calculate quantities like sample
# sizes, cluster counts, and time periods. These are used to auto-generate
# the cost and sample size functions.
#
# -----------------------------------------------------------------------------
# NAMING CONVENTION
# -----------------------------------------------------------------------------
#
#   {resource_type}_{stage}
#
# Where:
#   - resource_type: n (observations), clusters, periods, etc.
#   - stage: s1 (stage 1) or s2 (stage 2)
#
# Examples:
#   n_s1        : Observations in stage 1
#   n_s2        : Observations in stage 2
#   clusters_s1 : Clusters in stage 1
#   clusters_s2 : NEW clusters in stage 2
#   periods_s2  : Additional periods in stage 2
#
# -----------------------------------------------------------------------------
# FORMULA SYNTAX
# -----------------------------------------------------------------------------
#
# Use R formula syntax with ~ on the left:
#
#   resources = list(
#     n_s1       = ~ n_arms * k1 * m1,
#     n_s2       = ~ n_arms * k1 * m2 + n_arms * k2 * m2,
#     clusters_s1 = ~ n_arms * k1,
#     clusters_s2 = ~ n_arms * k2
#   )
#
# Available variables in formulas:
#   - All stage 1 parameters (k1, m1, K, t1, etc.)
#   - All stage 2 parameters (k2, m2, t2, r2, etc.)
#   - n_arms (from design spec)
#   - rho (cost ratio, when evaluating costs)
#   - Any outputs from model_builder
#
# -----------------------------------------------------------------------------
# HOW RESOURCES ARE USED
# -----------------------------------------------------------------------------
#
# 1. COST FUNCTION (auto-generated):
#
#    cost = sum(resource_qty * weight)
#
#    Default weights:
#      - n (observations): weight = 1
#      - clusters: weight = rho
#      - periods: weight = rho
#
#    For lambda-calibration, only stage 2 resources are counted
#    For fixed comparison, all resources are counted
#
# 2. SAMPLE SIZE FUNCTION (auto-generated):
#
#    n_stage1 = eval(n_s1 formula)
#    n_stage2 = eval(n_s2 formula)  # Only if continuing
#    n_total = n_stage1 + n_stage2
#

# =============================================================================
# COST STRUCTURE SPECIFICATION
# =============================================================================
#
# The cost_structure defines how resources translate to costs:
#
#   cost_structure = list(
#     weights = c(n = 1, clusters = "rho", periods = "rho"),
#     stage2_resources = c("n_s2", "clusters_s2")
#   )
#
# - weights: Named vector mapping resource types to cost weights
#   - Numeric values used directly
#   - Character values looked up as variables (e.g., "rho")
#
# - stage2_resources: Which resources to include for lambda-calibration
#   - Only these are used when finding optimal lambda
#   - Total cost uses all resources
#

# =============================================================================
# CUSTOM COST FUNCTION (optional)
# =============================================================================
#
# If the auto-generated cost function doesn't fit, write your own:
#
# -----------------------------------------------------------------------------
# SIGNATURE
# -----------------------------------------------------------------------------
#
#   cost_fn <- function(model_output, rho) { ... }
#
# -----------------------------------------------------------------------------
# INPUTS
# -----------------------------------------------------------------------------
#
# model_output: The list returned by your model_builder
#   - Contains all the fields you included (n1, n2, k1, k2, etc.)
#   - Access as: model_output$n2, model_output$k2, etc.
#
# rho: The cost ratio (from design$rho)
#
# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------
#
# Return a single scalar: the STAGE 2 cost (not total cost)
#
# This is used in the criterion: CP(z1) - lambda * cost
# The lambda-calibration finds the lambda that achieves target power.
#
# -----------------------------------------------------------------------------
# EXAMPLE
# -----------------------------------------------------------------------------
#
#   my_cost_fn <- function(model_output, rho) {
#     # Observations cost 1 each
#     n_cost <- model_output$n2
#
#     # Clusters cost rho each
#     k_cost <- 2 * model_output$k2  # 2 arms
#
#     # Return total stage 2 cost
#     n_cost + rho * k_cost
#   }
#

# =============================================================================
# CUSTOM SAMPLE SIZE FUNCTION (optional)
# =============================================================================
#
# If the auto-generated sample size function doesn't fit, write your own:
#
# -----------------------------------------------------------------------------
# SIGNATURE
# -----------------------------------------------------------------------------
#
#   sample_size_fn <- function(opt_designs, results, rho = 30) { ... }
#
# -----------------------------------------------------------------------------
# INPUTS
# -----------------------------------------------------------------------------
#
# opt_designs: Data.frame with one row per quadrature point (z1 value)
#   Columns include:
#     - z1: Stage 1 test statistic value
#     - continue: TRUE if trial continues at this z1
#     - efficacy_stop: TRUE if stopping for efficacy
#     - futility_stop: TRUE if stopping for futility
#     - cp: Conditional power at this z1
#     - All stage 2 design parameters (m2, k2, etc.)
#     - NA for design params when not continuing
#
# results: Full results object from find_lambda_for_power
#   Key paths:
#     - results$models$list[[1]]: First model output (has stage 1 params)
#     - results$models$design_grid: Stage 2 design grid
#     - results$quadrature$weights: Quadrature weights for integration
#
# rho: Cost ratio (provide default)
#
# -----------------------------------------------------------------------------
# REQUIRED OUTPUT
# -----------------------------------------------------------------------------
#
# Return a list with:
#
#   n_stage1    : Scalar - observations in stage 1
#   n_stage2    : Vector (length = nrow(opt_designs)) - observations per z1
#   n_total_max : Scalar - maximum possible total observations
#   n_total_min : Scalar - minimum total (usually = n_stage1)
#
# -----------------------------------------------------------------------------
# OPTIONAL OUTPUT
# -----------------------------------------------------------------------------
#
#   metrics: Named list of additional summary quantities
#     - Added to the summary table
#     - Examples: E_clusters, E_cost, max_cost, etc.
#
# -----------------------------------------------------------------------------
# EXAMPLE
# -----------------------------------------------------------------------------
#
#   my_sample_size_fn <- function(opt_designs, results, rho = 30) {
#
#     # Get stage 1 parameters
#     model1 <- results$models$list[[1]]
#     k1 <- model1$k1
#     m1 <- model1$m1
#     n_arms <- 2
#
#     # Stage 1 sample size (fixed)
#     n_stage1 <- n_arms * k1 * m1
#
#     # Stage 2 sample sizes (vector, one per z1)
#     m2 <- ifelse(is.na(opt_designs$m2), 0, opt_designs$m2)
#     k2 <- ifelse(is.na(opt_designs$k2), 0, opt_designs$k2)
#
#     n_stage2 <- ifelse(
#       opt_designs$continue,
#       n_arms * k1 * m2 + n_arms * k2 * m2,
#       0
#     )
#
#     # Maximum from design grid
#     grid <- results$models$design_grid
#     max_m2 <- max(grid$m2)
#     max_k2 <- max(grid$k2)
#     n_stage2_max <- n_arms * k1 * max_m2 + n_arms * max_k2 * max_m2
#
#     # Cost calculations
#     cost_stage1 <- n_stage1 + rho * (n_arms * k1)
#     cost_stage2 <- ifelse(
#       opt_designs$continue,
#       (n_arms * k1 * m2 + n_arms * k2 * m2) + rho * (n_arms * k2),
#       0
#     )
#
#     # Expected cost (using quadrature weights)
#     weights <- results$quadrature$weights
#     E_cost <- cost_stage1 + sum(cost_stage2 * weights)
#     max_cost <- cost_stage1 + n_stage2_max + rho * (n_arms * max_k2)
#
#     list(
#       n_stage1 = n_stage1,
#       n_stage2 = n_stage2,
#       n_total_max = n_stage1 + n_stage2_max,
#       n_total_min = n_stage1,
#       metrics = list(
#         E_cost = E_cost,
#         max_cost = max_cost,
#         E_clusters = n_arms * k1 + sum(n_arms * k2 * weights * opt_designs$continue)
#       )
#     )
#   }
#

# =============================================================================
# COMPLETE WORKED EXAMPLE: STAGGERED ENROLLMENT DESIGN
# =============================================================================
#
# A design where enrollment is staggered over time, with the option to
# extend enrollment in stage 2.

staggered_enrollment_example <- function() {

  cat("
================================================================================
EXAMPLE: Staggered Enrollment Design
================================================================================

Design:
  - K clusters total (fixed)
  - Stage 1: Enroll clusters over t1 periods at rate r1 per period
  - Stage 2: Optionally extend enrollment by t2 periods at rate r2
  - Outcome measured at fixed time after enrollment

Parameters:
  - Stage 1 (fixed): K, t1, r1, m

  - Stage 2 (adaptive): t2 (additional periods), r2 (enrollment rate)

")

  # -------------------------------------------------------------------------
  # Step 1: Write the model builder
  # -------------------------------------------------------------------------

  staggered_model_builder <- function(design_params, fixed_params) {

    # Extract parameters
    K <- fixed_params$K
    t1 <- fixed_params$t1
    r1 <- fixed_params$r1
    m <- fixed_params$m
    t2 <- design_params$t2
    r2 <- design_params$r2
    icc <- fixed_params$icc
    cac <- fixed_params$cac %||% 0.8
    delta <- fixed_params$delta
    cache <- fixed_params$cache

    T_total <- t1 + t2

    # Calculate enrollment times
    enroll_time <- numeric(K)
    enrolled <- 0

    for (t in 1:t1) {
      n_new <- min(r1, K - enrolled)
      if (n_new > 0) {
        enroll_time[(enrolled + 1):(enrolled + n_new)] <- t
        enrolled <- enrolled + n_new
      }
    }

    if (t2 > 0 && enrolled < K) {
      for (t in (t1 + 1):T_total) {
        n_new <- min(r2, K - enrolled)
        if (n_new > 0) {
          enroll_time[(enrolled + 1):(enrolled + n_new)] <- t
          enrolled <- enrolled + n_new
        }
        if (enrolled >= K) break
      }
    }

    # Clusters not enrolled get time > T_total (won't be observed)
    enroll_time[enroll_time == 0] <- T_total + 1

    # Build data: each cluster observed once at enrollment time
    # Treatment is randomized at enrollment
    df <- data.frame(
      cl = 1:K,
      t = enroll_time,
      trt = rep(0:1, length.out = K)  # Alternating assignment
    )
    df <- df[df$t <= T_total, ]  # Only enrolled clusters
    df$n <- m

    # Stage indices
    idx1 <- which(df$t <= t1)
    idx2 <- which(df$t > t1)

    n1 <- length(idx1)
    n2 <- length(idx2)

    # Build model
    v1 <- icc
    cache_key <- paste(K, T_total, r2, sep = "_")

    if (is.null(cache$mod) || cache$cache_key != cache_key) {
      cache$mod <- Model$new(
        ~ trt + (1|gr(cl)),
        data = df,
        family = gaussian(),
        mean = c(0, delta),
        covariance = v1,
        weights = df$n,
        var_par = 1 - icc
      )
      cache$cache_key <- cache_key
    } else {
      glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    }

    # Compute V
    X <- cache$mod$mean$X
    Z <- cache$mod$covariance$Z
    D <- cache$mod$covariance$D

    Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
    D_sp <- Matrix::Matrix(D, sparse = TRUE)
    S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) +
      Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
    V <- as.matrix(S)

    # Efficient score decomposition
    eff <- efficient_score_decomposition(X, V, idx1, idx2, j = 2)

    list(
      # Required
      I1_eff = eff$I1_eff,
      I2_eff = eff$I2_eff,
      I_eff = eff$I_eff,
      w1 = sqrt(eff$I1_eff / eff$I_eff),
      w2 = sqrt(eff$I2_eff / eff$I_eff),
      b1 = delta,
      # Resource tracking
      n1 = n1 * m,
      n2 = n2 * m,
      clusters_s1 = n1,
      clusters_s2 = n2,
      K = K, t1 = t1, t2 = t2, r1 = r1, r2 = r2, m = m,
      T_total = T_total
    )
  }

  # -------------------------------------------------------------------------
  # Step 2: Create design specification
  # -------------------------------------------------------------------------

  staggered_spec <- crt_design_spec(
    stage1_params = c("K", "t1", "r1", "m"),
    stage2_params = c("t2", "r2"),
    resources = list(
      n_s1 = ~ clusters_s1 * m,  # Uses model output
      n_s2 = ~ clusters_s2 * m,
      clusters_s1 = ~ clusters_s1,
      clusters_s2 = ~ clusters_s2,
      periods_s1 = ~ t1,
      periods_s2 = ~ t2
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho", periods = "rho"),
      stage2_resources = c("n_s2", "clusters_s2", "periods_s2")
    ),
    model_builder = staggered_model_builder,
    n_arms = 1,  # Not really arms in this design
    design_type = "staggered_enrollment"
  )

  # -------------------------------------------------------------------------
  # Step 3: Create design object
  # -------------------------------------------------------------------------

  staggered_design <- structure(
    list(
      spec = staggered_spec,
      fixed_params = list(icc = 0.05, cac = 0.8, delta = 0.3),
      stage1_grid = data.frame(K = 30, t1 = 4, r1 = 4, m = 25),
      stage2_grid_fn = function(s1) {
        expand.grid(t2 = 0:6, r2 = c(2, 4, 6, 8))
      },
      rho = 50  # Cost per additional period/cluster
    ),
    class = c("staggered_crt", "crt_design")
  )

  cat("
Design object created. To analyse:

  results <- adaptive_analysis(staggered_design, target_power = 0.8)
  summary(results)
  plot(results, type = 'decision')

")

  invisible(staggered_design)
}


# =============================================================================
# HELP FUNCTION
# =============================================================================

#' Print help on custom design creation
#'
#' @param topic Optional topic: "model_builder", "resources", "cost_fn",
#'        "sample_size_fn", "example", or NULL for overview
#' @export
custom_design_help <- function(topic = NULL) {

  if (is.null(topic)) {
    cat("
================================================================================
                     CUSTOM DESIGN GUIDE FOR adaptiveCRT
================================================================================

To create a custom adaptive trial design, you need:

  1. model_builder  - Function that computes efficient information (REQUIRED)
  2. resources      - Formulas for sample sizes, clusters, etc. (REQUIRED)
  3. cost_fn        - Custom cost function (optional, auto-generated by default
)
  4. sample_size_fn - Custom sample size function (optional, auto-generated)

QUICK START:
------------
  # See complete working example:
  custom_design_help('example')

DETAILED TOPICS:
----------------
  custom_design_help('model_builder')   # How to write the model function

  custom_design_help('resources')       # How to specify resource formulas
  custom_design_help('cost_fn')         # How to write custom cost function
  custom_design_help('sample_size_fn')  # How to write custom sample size function
  custom_design_help('example')         # Complete worked example

BASIC TEMPLATE:
---------------
  # 1. Write model builder
  my_model <- function(design_params, fixed_params) {
    # ... compute efficient information ...
    list(I1_eff=..., I2_eff=..., I_eff=..., w1=..., w2=..., b1=...)
  }

  # 2. Create specification
  my_spec <- crt_design_spec(
    stage1_params = c('k1', 'm1'),
    stage2_params = c('k2', 'm2'),
    resources = list(
      n_s1 = ~ 2 * k1 * m1,
      n_s2 = ~ 2 * k1 * m2 + 2 * k2 * m2,
      clusters_s1 = ~ 2 * k1,
      clusters_s2 = ~ 2 * k2
    ),
    model_builder = my_model
  )

  # 3. Create design object
  my_design <- structure(list(
    spec = my_spec,
    fixed_params = list(icc = 0.05, delta = 0.25),
    stage1_grid = expand.grid(k1 = 8:12, m1 = c(20, 30)),
    stage2_grid_fn = function(s1) expand.grid(m2 = 10:50, k2 = 0:4),
    rho = 30
  ), class = c('my_crt', 'crt_design'))

  # 4. Analyse
  results <- adaptive_analysis(my_design, target_power = 0.8)

================================================================================
")
  } else if (topic == "model_builder") {
    cat("
================================================================================
                            MODEL_BUILDER FUNCTION
================================================================================

SIGNATURE:
----------
  model_builder <- function(design_params, fixed_params) { ... }

INPUTS:
-------
  design_params: Single-row data.frame with stage 2 parameters
    - Example: data.frame(m2 = 30, k2 = 2)
    - Access as: design_params$m2

  fixed_params: List containing:
    - Stage 1 parameters (k1, m1, etc. - merged from stage1_grid)
    - Constants (icc, cac, delta)
    - cache: Environment for model caching (IMPORTANT!)
    - Access as: fixed_params$k1, fixed_params$icc, fixed_params$cache

REQUIRED OUTPUTS:
-----------------
  Return a list with:

    I1_eff  : Stage 1 efficient information (scalar)
    I2_eff  : Stage 2 efficient information (scalar)
    I_eff   : Total efficient information (scalar)
    w1      : sqrt(I1_eff / I_eff)
    w2      : sqrt(I2_eff / I_eff)
    b1      : Treatment effect (= delta)

OPTIONAL OUTPUTS (for resource tracking):
-----------------------------------------
    n1, n2        : Observations per stage
    k1, k2, m1, m2: Design parameters
    df_s1, df_full: Degrees of freedom for t-distribution
    ... any quantities needed by cost/sample size functions

COMPUTING EFFICIENT INFORMATION:
--------------------------------
  # Use the helper function:
  eff <- efficient_score_decomposition(X, V, idx1, idx2, j = j_treatment)

  # Where:
  #   X = fixed effects design matrix
  #   V = marginal covariance matrix
  #   idx1 = row indices for stage 1 data
  #   idx2 = row indices for stage 2 data
  #   j = column index of treatment in X

CACHING (CRITICAL FOR SPEED):
-----------------------------
  cache <- fixed_params$cache
  cache_key <- paste(k1, k2, sep = '_')

  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- Model$new(...)
    cache$cache_key <- cache_key
  } else {
    # Just update weights/parameters
    glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
  }

================================================================================
")
  } else if (topic == "resources") {
    cat("
================================================================================
                           RESOURCES SPECIFICATION
================================================================================

Resources are formulas that define quantities like sample sizes and cluster
counts. They're used to auto-generate cost and sample size functions.

NAMING CONVENTION:
------------------
  {type}_{stage}

  Examples:
    n_s1        : Observations in stage 1
    n_s2        : Observations in stage 2
    clusters_s1 : Clusters in stage 1
    clusters_s2 : NEW clusters in stage 2
    periods_s2  : Additional time periods in stage 2

SYNTAX:
-------
  Use R formula notation with ~ on the left:

  resources = list(
    n_s1        = ~ n_arms * k1 * m1,
    n_s2        = ~ n_arms * k1 * m2 + n_arms * k2 * m2,
    clusters_s1 = ~ n_arms * k1,
    clusters_s2 = ~ n_arms * k2
  )

AVAILABLE VARIABLES:
--------------------
  - All stage 1 parameters (k1, m1, K, t1, etc.)
  - All stage 2 parameters (k2, m2, t2, r2, etc.)
  - n_arms (from design spec)
  - rho (when computing costs)
  - Any fields from model_builder output

COST STRUCTURE:
---------------
  cost_structure = list(
    weights = c(n = 1, clusters = 'rho', periods = 'rho'),
    stage2_resources = c('n_s2', 'clusters_s2')
  )

  - weights: Maps resource types to cost multipliers
    - Numeric: used directly
    - Character: looked up as variable (e.g., 'rho')

  - stage2_resources: Which resources count for lambda-calibration
    - Only these are used in the CP - lambda*cost criterion

================================================================================
")
  } else if (topic == "cost_fn") {
    cat("
================================================================================
                          CUSTOM COST FUNCTION
================================================================================

Only needed if auto-generated cost doesn't fit your design.

SIGNATURE:
----------
  cost_fn <- function(model_output, rho) { ... }

INPUTS:
-------
  model_output: List returned by your model_builder
    - Contains all fields you included (n1, n2, k1, k2, etc.)
    - Access as: model_output$n2

  rho: Cost ratio from design$rho

OUTPUT:
-------
  Return scalar: STAGE 2 cost only (not total)

  Used in criterion: CP(z1) - lambda * cost

EXAMPLE:
--------
  my_cost_fn <- function(model_output, rho) {
    n_cost <- model_output$n2
    k_cost <- 2 * model_output$k2
    n_cost + rho * k_cost
  }

================================================================================
")
  } else if (topic == "sample_size_fn") {
    cat("
================================================================================
                       CUSTOM SAMPLE SIZE FUNCTION
================================================================================

Only needed if auto-generated sample size doesn't fit your design.

SIGNATURE:
----------
  sample_size_fn <- function(opt_designs, results, rho = 30) { ... }

INPUTS:
-------
  opt_designs: Data.frame with one row per z1 quadrature point
    Columns:
      - z1: Stage 1 test statistic
      - continue: TRUE if continuing to stage 2
      - efficacy_stop, futility_stop: Early stopping indicators
      - All stage 2 parameters (NA when not continuing)

  results: Full results object
    Key paths:
      - results$models$list[[1]]: First model (has stage 1 params)
      - results$models$design_grid: Stage 2 grid
      - results$quadrature$weights: For computing expectations

REQUIRED OUTPUT:
----------------
  list(
    n_stage1    = scalar,
    n_stage2    = vector (length = nrow(opt_designs)),
    n_total_max = scalar,
    n_total_min = scalar
  )

OPTIONAL OUTPUT:
----------------
  metrics = list(
    E_cost = ...,
    max_cost = ...,
    E_clusters = ...
  )

  These appear in the summary table.

================================================================================
")
  } else if (topic == "example") {
    staggered_enrollment_example()
  } else {
    cat("Unknown topic:", topic, "\n")
    cat("Available topics: model_builder, resources, cost_fn, sample_size_fn, example\n")
  }

  invisible(NULL)
}
