# =============================================================================
# IMPALA-style adaptive crossover trial: from design to analysis
# =============================================================================
# A two-sequence cluster crossover trial with a binary outcome. Stage 1 is the
# first full crossover (two periods). At the interim the cluster-period size for
# the second crossover is chosen. The script designs the trial, simulates one
# realisation of it and analyses that realisation with the package functions.
#
# The model builder and design constructor are specific to this design, so they
# are defined here. Everything else (nuisance estimation, score statistics,
# interim decision, effect estimates) comes from the package.
#
# Run it with:
#   source(system.file("examples", "crossover_design_to_analysis.R", package = "acrt"))

library(acrt)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
set.seed(20260902)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
FAMILY     <- "binomial"
BASELINE_P <- 0.15               # planned control-arm risk
DELTA      <- -0.03              # planned risk difference

K1     <- 15L                    # stage 1 clusters per sequence
M11    <- 75L                    # cluster size, period 1
M12    <- 75L                    # cluster size, period 2
T2     <- 2L                     # stage 2 periods: the second full crossover
K2MAX  <- 0L                     # maximum new clusters per sequence at stage 2
M2GRID <- seq(30L, 80L, by = 10L)
T_MAX  <- 2L + T2

RHO   <- 300                     # cost per cluster
RHO_T <- 30                      # cost per period

ICC_P <- 0.10                    # planning ICC
CAC_P <- 0.80                    # planning lag-1 cluster autocorrelation

MARGINAL_TARGET <- TRUE          # BASELINE_P and DELTA are marginal quantities

# -----------------------------------------------------------------------------
# 1. Model builder and design constructor for the crossover design
# -----------------------------------------------------------------------------

b1_link_at <- function(p0, rd = DELTA) qlogis(p0 + rd) - qlogis(p0)

marginal_p <- function(eta, A, nodes = 32L) {
  gh <- statmod::gauss.quad(nodes, kind = "hermite")
  x <- gh$nodes * sqrt(2 * A); w <- gh$weights / sqrt(pi)
  sum(w * plogis(eta + x))
}

calibrate_link <- function(p0, rd, A) {
  if (!MARGINAL_TARGET || !is.finite(A) || A <= 0)
    return(list(eta0 = qlogis(p0), b1 = if (rd == 0) 0 else b1_link_at(p0, rd)))
  e0 <- uniroot(function(e) marginal_p(e, A) - p0,
                qlogis(p0) + c(-4, 4))$root
  if (rd == 0) return(list(eta0 = e0, b1 = 0))
  e1 <- uniroot(function(e) marginal_p(e, A) - (p0 + rd),
                qlogis(p0 + rd) + c(-4, 4))$root
  list(eta0 = e0, b1 = e1 - e0)
}

convert_to_glmm_scale <- function(family = c("gaussian", "binomial", "poisson"),
                                  baseline = 0,
                                  delta,
                                  icc,
                                  cac = NULL) {

  family <- match.arg(family)

  if (family == "gaussian") {

    beta <- c(baseline, delta)

    cov_pars <- if (!is.null(cac)) {
      c(icc * cac, icc * (1 - cac))
    } else {
      icc
    }

    return(list(
      beta = beta,
      cov_pars = cov_pars,
      var_par = 1 - icc,
      family = gaussian(),
      b1_original = delta,
      b1_link = delta,
      link = "identity"
    ))

  } else if (family == "binomial") {

    p0 <- baseline
    p1 <- baseline + delta

    if (p0 <= 0 || p0 >= 1)
      stop("baseline must lie in (0,1).")

    if (p1 <= 0 || p1 >= 1)
      stop("baseline + delta must lie in (0,1).")

    ## Approximate marginal-ICC -> latent variance mapping
    A <- (icc / (1 - icc)) / (p0 * (1 - p0))

    ## Exact marginal mean calibration conditional on A
    lk <- calibrate_link(p0, delta, A)

    beta <- c(lk$eta0, lk$b1)

    cov_pars <- if (!is.null(cac)) {
      c(A * cac, A * (1 - cac))
    } else {
      A
    }

    return(list(
      beta = beta,
      cov_pars = cov_pars,
      tau_sq = A,
      var_par = 1,
      family = binomial(),
      b1_original = delta,
      b1_link = lk$b1,
      link = "logit"
    ))

  } else if (family == "poisson") {

    b0 <- log(baseline)
    rate1 <- baseline + delta

    if (rate1 <= 0)
      stop("baseline + delta must be positive.")

    b1 <- log(rate1) - b0
    tau_sq <- (icc / (1 - icc)) / baseline

    cov_pars <- if (!is.null(cac)) {
      c(tau_sq * cac, tau_sq * (1 - cac))
    } else {
      tau_sq
    }

    return(list(
      beta = c(b0, b1),
      cov_pars = cov_pars,
      tau_sq = tau_sq,
      var_par = 1,
      family = poisson(),
      b1_original = delta,
      b1_link = b1,
      link = "log"
    ))
  }
}

crossover_model_builder <- function(design_params, fixed_params) {
  cache <- fixed_params$cache
  k1 <- fixed_params$k1
  m11 <- fixed_params$m11
  m12 <- fixed_params$m12
  k2 <- design_params$k2
  m2 <- design_params$m2
  m2b <- design_params$m2b %||% design_params$m2
  icc <- fixed_params$icc
  cac <- fixed_params$cac
  delta <- fixed_params$delta
  cache <- fixed_params$cache
  t2 <- design_params$t2
  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0

  # Convert parameters to appropriate scale
  glm_params <- convert_to_glmm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = cac
  )
  # Stage 1: crossover (2 periods)
  df1 <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k1, ")) * t(2)")))
  df1$int <- df1$int - 1
  df1$t <- df1$t - 1
  # Treatment: arm 1 gets trt in period 2, arm 2 gets trt in period 1
  df1$trt <- ifelse((df1$int == 0 & df1$t == 1) | (df1$int == 1 & df1$t == 0), 1, 0)
  df1$n <- m11
  df1[df1$t == 1, "n"] <- m12
  if(t2 == 1){
    # Stage 2 period (period 3)
    df2_base <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
    df2_base$int <- df2_base$int - 1
    df2_base$t <- 2
    df2_base$trt <- df2_base$int  # Same as their original arm
    df2_base$n <- m2
  } else {
    df2_base <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k1, ")) * t(2)")))
    df2_base$int <- df2_base$int - 1
    df2_base$trt <- ifelse((df2_base$int == 0 & df2_base$t == 2) | (df2_base$int == 1 & df2_base$t == 1), 1, 0)
    df2_base$n <- m2
    df2_base$t <- df2_base$t + 1
  }

  if (k2 > 0) {
    # New clusters get both periods
    df2_new <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k2, ")) * t(2)")))
    df2_new$int <- df2_new$int - 1
    df2_new$trt <- ifelse((df2_new$int == 0 & df2_new$t == 2) | (df2_new$int == 1 & df2_new$t == 1), 1, 0)
    df2_new$n <- m2
    df2_new$t <- df2_new$t + 1
    df2_new$cl <- df2_new$cl + max(df1$cl)
    df <- rbind(df1, df2_base, df2_new)
  } else {
    df <- rbind(df1, df2_base)
  }
  # Number of time periods
  n_periods <- length(unique(df$t))

  # Mean parameters: intercept + (n_periods - 1) time effects + treatment
  beta_full <- c(glm_params$beta[1],           # Intercept
                 rep(0, n_periods - 1),         # Time effects
                 glm_params$beta[2]         )   # Treatment

  # gr(cl)*ar0(t) takes (variance, AR parameter). convert_to_glm_scale()
  # returns the nested-exchangeable split (tau2*cac, tau2*(1-cac)), so
  # recombine it once and use it for every family and on cache hits.
  cov_ar <- c(sum(glm_params$cov_pars), cac)

  # Cache key includes family
  cache_key <- paste(k1, k2, t2, family, sep = "_")
  if (is.null(cache$mod) || cache$cache_key != cache_key) {

    if (family == "gaussian") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = cov_ar,
        weights = df$n,
        var_par = glm_params$var_par
      )
    } else if (family == "binomial") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = cov_ar,
        trials = df$n
      )
    } else if (family == "poisson") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = cov_ar,
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
    cache$mod$update_parameters(cov.pars = cov_ar, mean.pars = beta_full)
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
    j_trt <- ncol(X)
  }
  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)
  # Degrees of freedom
  n_cluster_periods_s1 <- nrow(df1)
  n_cluster_periods_full <- nrow(df)
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
    m11 = m11,
    m12 = m12,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2,
    t2 = t2,
    df_s1 = df_s1,
    df_s2 = (nrow(df) - nrow(df1)) - (n_periods - 2) - 1,
    df_full = df_full,
    family = family,
    link = glm_params$link
  )
}

make_crossover_model_fn <- function() {
  cache <- new.env()
  cache$mod <- NULL
  cache$cache_key <- NULL

  function(design_params, fixed_params) {
    fixed_params$cache <- cache
    crossover_model_builder(design_params, fixed_params)
  }
}

crossover_crt <- function(icc, cac = 0.8, delta,
                          k1, m11, m12,
                          k2 = 0:4, m2 = NULL, t2 = 1:2,
                          rho = 30, rho_t = 10,
                          family = "gaussian", baseline = NULL) {

  spec <- crt_design_spec(
    stage1_params = c("k1", "m11", "m12"),
    stage2_params = c("k2", "m2", "t2"),
    resources = list(
      n_s1 = ~ n_arms * k1 * (m11 + m12),        # 2 periods in stage 1
      n_s2 = ~ n_arms * k1 * m2 * t2 + n_arms * k2 * m2 * 2,  # new clusters run two periods
      clusters_s1 = ~ n_arms * k1,
      clusters_s2 = ~ n_arms * k2,
      t_s1 = ~ 2,
      t_s2 = ~ max(t2, 2 * (k2 > 0))
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho", t = "rho_t"),
      stage2_resources = c("n_s2", "clusters_s2", "t_s2")
    ),
    model_builder = crossover_model_builder,
    n_arms = 2,
    design_type = "crossover"
  )
  structure(
    list(
      spec = spec,
      fixed_params = list(
        icc = icc,
        cac = cac,
        delta = delta,
        family = family,
        baseline = baseline %||% 0,
        rho_t = rho_t
      ),
      stage1_grid = expand.grid(k1 = k1, m11 = m11, m12 = m12),
      stage2_grid_fn = function(s1) expand.grid(t2 = t2, k2 = k2, m2 = m2),
      rho = rho
    ),
    class = c("crossover_crt", "crt_design")
  )
}

# -----------------------------------------------------------------------------
# 2. Design the trial
# -----------------------------------------------------------------------------

design <- crossover_crt(
  icc      = ICC_P,
  cac      = CAC_P,
  delta    = DELTA,
  baseline = BASELINE_P,
  k1       = K1,
  m11      = M11,
  m12      = M12,
  k2       = 0:K2MAX,
  m2       = M2GRID,
  t2       = T2,
  rho      = RHO,
  rho_t    = RHO_T,
  family   = FAMILY
)
design$spec$model_builder <- make_crossover_model_fn()

planned <- adaptive_analysis(design, target_power = 0.8, tol = 0.005,
                             method = "lambda", verbose = FALSE)
print(planned)

# -----------------------------------------------------------------------------
# 3. Fix the analysis plan
# -----------------------------------------------------------------------------
# The weights, boundaries, family and planned CAC come from the design. The
# schedule follows the treatment rule in the model builder. Sequences 3 and 4
# are for any clusters added at stage 2. The design is balanced, so no period
# effects are included. If randomisation were stratified, add strata = "country".

schedule <- rbind(
  c( 0,  1, 0, 1),     # sequence 1: control first
  c( 1,  0, 1, 0),     # sequence 2: intervention first
  c(NA, NA, 0, 1),     # sequence 3: added at stage 2
  c(NA, NA, 1, 0)      # sequence 4: added at stage 2
)

plan <- analysis_plan(planned, schedule = schedule, stage1 = 1:2)
print(plan)

# -----------------------------------------------------------------------------
# 4. Simulate stage 1
# -----------------------------------------------------------------------------
# Cluster random effects follow an AR(1) process over the periods, with the
# intercept and treatment effect calibrated to the planned marginal risks.

truth <- convert_to_glmm_scale(FAMILY, BASELINE_P, DELTA, ICC_P, CAC_P)
A <- truth$tau_sq
n_clusters <- 2L * (K1 + K2MAX)
cluster_effects <- matrix(rnorm(n_clusters * T_MAX), n_clusters, T_MAX) %*%
  chol(A * CAC_P^abs(outer(seq_len(T_MAX), seq_len(T_MAX), "-")))

simulate_periods <- function(alloc, period, m) {
  on_intervention <- schedule[alloc$sequence, period]
  logit <- truth$beta[1] + truth$b1_link * on_intervention +
    cluster_effects[alloc$cluster, period]
  data.frame(cluster = alloc$cluster, period = period, m = m,
             y = rbinom(nrow(alloc), m, plogis(logit)))
}

alloc <- data.frame(cluster = seq_len(2L * K1),
                    sequence = rep(1:2, each = K1))

data_stage1 <- rbind(simulate_periods(alloc, 1L, M11),
                     simulate_periods(alloc, 2L, M12))

# -----------------------------------------------------------------------------
# 5. Interim analysis
# -----------------------------------------------------------------------------
# Futility is non-binding, so override_stop = TRUE still returns the stage 2
# design when the futility rule fires. Only an efficacy stop ends the trial here.

interim <- interim_statistic(plan, data_stage1, alloc)
print(interim)

decision <- interim_decision(interim, override_stop = TRUE, verbose = FALSE)
print(decision)

# -----------------------------------------------------------------------------
# 6. Stage 2 and final analysis, or stop at the interim
# -----------------------------------------------------------------------------

if (decision$decision != "efficacy_stop") {

  k2 <- decision$stage2_design$k2
  m2 <- decision$stage2_design$m2

  alloc_final <- alloc
  if (k2 > 0) {
    alloc_final <- rbind(
      alloc,
      data.frame(cluster = 2L * K1 + seq_len(2L * k2),
                 sequence = rep(3:4, each = k2))
    )
  }

  data_final <- rbind(data_stage1,
                      simulate_periods(alloc_final, 3L, m2),
                      simulate_periods(alloc_final, 4L, m2))

  final <- final_statistic(interim, data_final, alloc = alloc_final)
  print(final)

  estimates <- effect_estimates(final)

} else {

  estimates <- effect_estimates(interim)
}

# -----------------------------------------------------------------------------
# 7. Results
# -----------------------------------------------------------------------------
# The `adaptive` row is the primary result. Its interval is the repeated
# confidence interval. The other rows are descriptive.

print(estimates)
cat(sprintf("\nTrue log odds ratio (cluster-specific): %.3f\n", truth$b1_link))
cat(sprintf("True marginal risk difference: %.3f\n", DELTA))
