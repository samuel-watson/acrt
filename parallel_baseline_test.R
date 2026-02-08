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
parallel_baseline_crt <- function(icc, delta, k1, m1, m1b,
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
    stage1_params = c("k1", "m1", "m1b"),
    stage2_params = c("k2", "m2"),
    resources = list(
      n_s1 = ~ n_arms * k1 * (m1 + m1b),
      n_s2 = ~ n_arms * k1 * m2 + n_arms * k2 * m2,
      clusters_s1 = ~ n_arms * k1,
      clusters_s2 = ~ n_arms * k2
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho"),
      stage2_resources = c("n_s2", "clusters_s2")
    ),
    model_builder = parallel_baseline_model_builder,
    n_arms = n_arms,
    design_type = "parallel_baseline"
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
      stage1_grid = expand.grid(k1 = k1, m1 = m1, m1b = m1b),
      stage2_grid_fn = stage2_grid_fn,
      rho = rho
    ),
    class = c("parallel_baseline_crt", "crt_design")
  )
}

#' Model builder for parallel CRT (internal)
#' @keywords internal
parallel_baseline_model_builder <- function(design_params, fixed_params) {

  k1 <- fixed_params$k1
  m1 <- fixed_params$m1
  m1b <- fixed_params$m1b
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
  df2$t <- ifelse(m1b > 0, 2, 1)

  #baseline
  if(m1b > 0){
    df1b <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
    df1b$int <- 0
    df1b$n <- m1b
    df1b$t <- 0
    df1$t <- 1
    df1 <- rbind(df1, df1b)
  }

  # Stage 2: new clusters (if any)
  if (k2 > 0) {
    df2b <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k2, ")")))
    df2b$int <- df2b$int - 1
    df2b$cl <- df2b$cl + k1 * 2
    df2b$n <- m2b
    df2b$t <- ifelse(m1b > 0, 2, 1)
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
        ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, rep(0, ifelse(m1b > 0, 2, 1))),
        covariance = glm_params$cov_pars,
        weights = df$n,
        var_par = glm_params$var_par
      )
    } else if (family == "binomial") {
      cache$mod <- glmmrBase::Model$new(
        ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, rep(0, ifelse(m1b > 0, 2, 1))),
        covariance = glm_params$cov_pars,
        trials = df$n
      )
    } else if (family == "poisson") {
      # For Poisson, n represents exposure (e.g., person-time)
      cache$mod <- glmmrBase::Model$new(
        ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
        data = df,
        family = glm_params$family,
        mean = c(glm_params$beta, rep(0, ifelse(m1b > 0, 2, 1))),
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
    m1b = m1b,
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


# Also need factory function for caching
make_parallel_baseline_model_fn <- function() {
  cache <- new.env()
  cache$mod <- NULL
  cache$cache_key <- NULL

  function(design_params, fixed_params) {
    fixed_params$cache <- cache
    parallel_baseline_model_builder(design_params, fixed_params)
  }
}

# Define the design in ONE call - no custom functions needed!
design <- parallel_baseline_crt(
  icc = 0.018,           # Intra-cluster correlation
  cac = 0.97,            # Cluster autocorrelation
  delta = -0.005,         # Treatment effect
  baseline = 0.02,
  k1 = seq(20,32,by=4),    # Stage 1 clusters per arm (explore these)
  m1 = seq(100,1300,by=300),
  m1b = seq(0,1200,by=300),
  rho = 200,              # Cluster-to-individual cost ratio
  family = "binomial"
)

# Run the analysis - ONE function call
results <- adaptive_analysis(design, target_power = 0.9, verbose = TRUE, tol = 0.01)

# View summary
summary(results)

# Plot results
p1 <- plot(results, type = "EN")      # Expected sample size
p2 <- plot(results, type = "Nmax")    # Maximum sample size
p3 <- plot(results, type = "pareto", objectives = list(E_cost = "min", max_cost = "min"))  # Pareto frontier

p3 +
  geom_vline(xintercept = 201600, lty = 2) +
  geom_hline(yintercept = 201600, lty = 2)

# Define the design in ONE call - no custom functions needed!
design_single <- parallel_baseline_crt(
  icc = 0.018,           # Intra-cluster correlation
  cac = 0.97,            # Cluster autocorrelation
  delta = -0.005,         # Treatment effect
  baseline = 0.02,
  k1 = 32,    # Stage 1 clusters per arm (explore these)
  m1 = 400,
  m1b = 1200,
  rho = 200,              # Cluster-to-individual cost ratio
  family = "binomial"
)

# Analyse
results_single <- adaptive_analysis(design_single, target_power = 0.9, tol = 0.01)

# Summary shows more detail for single designs
summary(results_single)

# Decision rules
plot(results_single, type = "decision")
plot_decision_rules(results_single, design_var = "t2")
# Get decision rules as data frame
rules <- get_decision_rules(results_single)
print(rules[seq(1, nrow(rules), by = 5), c("z1", "continue", "cp", "m2", "k2")])


# Find Pareto-optimal designs (minimise E[N] and N_max)
pareto <- find_pareto(results, objectives = list(E_cost = "min", max_cost = "min"))
print(pareto[, c("k1", "m1", "E_N", "N_max", "power_stage1")])

# Different objectives: minimise cost, maximise early stopping
pareto2 <- find_pareto(results,
                       objectives = list(E_cost = "min", prob_efficacy = "max"))
print(pareto2[, c("k1", "m1", "E_cost", "prob_efficacy")])

# With constraints
pareto3 <- find_pareto(results,
                       objectives = list(E_N = "min", N_max = "min"),
                       constraints = list(power_stage1 = c(0.3, 0.6)))
print(pareto3[, c("k1", "m1", "E_N", "N_max", "power_stage1")])
