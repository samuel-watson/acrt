# =============================================================================
# Estimation and inference for completed or interim adaptive CRTs
# =============================================================================
# Stage-wise combination score test, point estimates and repeated confidence
# intervals, computed from trial data using the design object.
#
# Workflow:
#   plan <- analysis_plan(planned_results, schedule = ..., stage1 = ...)
#   int  <- interim_statistic(plan, data_interim, alloc)
#   dec  <- interim_decision(int)
#   fin  <- final_statistic(int, data_final)
#   est  <- effect_estimates(fin)
#
# Every statistic is evaluated at a hypothesised treatment effect beta0 using
# nuisance values locked at the interim, with stratum intercepts re-profiled
# under beta0. The analysis is deterministic: it uses no random numbers.

# =============================================================================
# Family definitions
# =============================================================================

.estimation_family <- function(family) {
  switch(
    family,
    binomial = list(
      name = "binomial",
      grid = TRUE,
      labels = c(link = "log odds ratio", ratio = "odds ratio",
                 difference = "risk difference"),
      has_ratio = TRUE,
      linkinv = stats::plogis,
      loglik = function(H, y, m) {
        y * stats::plogis(H, log.p = TRUE) +
          (m - y) * stats::plogis(-H, log.p = TRUE)
      },
      working = function(y, m, sigma2) {
        yc <- pmin(pmax(y / m, 0.5 / m), 1 - 0.5 / m)
        list(ell = stats::qlogis(yc), w = 1 / (m * yc * (1 - yc)))
      },
      info = function(y, m) {
        yc <- pmin(pmax(y / m, 0.5 / m), 1 - 0.5 / m)
        m * yc * (1 - yc)
      },
      clamp_mean = function(p, m, control) .clamp(p, control$p_bounds),
      icc_scale = function(mean, sigma2) mean * (1 - mean),
      marginal = function(eta, A, gh) {
        vapply(eta, function(e) {
          sum(gh$weights * stats::plogis(e + gh$nodes * sqrt(A)))
        }, numeric(1))
      },
      dmarginal = function(eta, A, gh) {
        drop(stats::dlogis(outer(eta, gh$nodes * sqrt(A), "+")) %*% gh$weights)
      },
      intercept = function(mean, share1, beta, A, gh) {
        marginal <- function(e) {
          sum(gh$weights * stats::plogis(e + gh$nodes * sqrt(A)))
        }
        stats::uniroot(
          function(e) {
            (1 - share1) * marginal(e) + share1 * marginal(e + beta) - mean
          },
          stats::qlogis(mean) - share1 * beta + c(-4, 4),
          extendInt = "yes", tol = 1e-10
        )$root
      }
    ),
    poisson = list(
      name = "poisson",
      grid = TRUE,
      labels = c(link = "log rate ratio", ratio = "rate ratio",
                 difference = "rate difference"),
      has_ratio = TRUE,
      linkinv = exp,
      loglik = function(H, y, m) y * H - m * exp(H),
      working = function(y, m, sigma2) {
        yc <- pmax(y, 0.5)
        list(ell = log(yc / m), w = 1 / yc)
      },
      info = function(y, m) pmax(y, 0.5),
      clamp_mean = function(p, m, control) max(p, 0.5 / m),
      icc_scale = function(mean, sigma2) mean,
      marginal = function(eta, A, gh) exp(eta + A / 2),
      dmarginal = function(eta, A, gh) exp(eta + A / 2),
      intercept = function(mean, share1, beta, A, gh) {
        log(mean) - A / 2 - log((1 - share1) + share1 * exp(beta))
      }
    ),
    gaussian = list(
      name = "gaussian",
      grid = FALSE,
      labels = c(link = "mean difference", ratio = NA_character_,
                 difference = "mean difference"),
      has_ratio = FALSE,
      linkinv = identity,
      loglik = NULL,
      working = function(y, m, sigma2) list(ell = y / m, w = sigma2 / m),
      clamp_mean = function(p, m, control) p,
      icc_scale = function(mean, sigma2) 1 / sigma2,
      marginal = function(eta, A, gh) eta,
      dmarginal = function(eta, A, gh) rep(1, length(eta)),
      intercept = function(mean, share1, beta, A, gh) mean - share1 * beta
    ),
    stop("family must be one of 'binomial', 'poisson' or 'gaussian'")
  )
}

.clamp <- function(x, b) min(max(x, b[1]), b[2])

.estimation_control <- function(control) {
  defaults <- list(
    cac_bounds = c(0, 0.99),
    icc_floor = 1e-4,
    p_bounds = c(0.01, 0.99),
    gh_nodes = 32L,
    span = 6,
    grid_min = 80L,
    grid_max = 2000L
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control setting(s): ", paste(unknown, collapse = ", "))
  }
  utils::modifyList(defaults, control)
}

# =============================================================================
# Analysis plan
# =============================================================================

#' Specify the analysis of an adaptive cluster randomised trial
#'
#' Collects everything fixed before the trial data are seen: the outcome
#' family, the treatment schedule, the combination weights and the boundaries.
#' When `planned` is supplied, these are taken from the design analysis so that
#' the trial is analysed with exactly the values it was designed with.
#'
#' @param planned An `adaptive_results` object for a single stage 1 design, from
#'   [adaptive_analysis()]. Supplies `family`, `w1`, `c1`, `c2` and `cac_plan`,
#'   and enables [interim_decision()]. If `NULL`, at least `family`, `schedule`,
#'   `stage1`, `w1` and `cac_plan` must be given.
#' @param schedule Numeric matrix with one row per sequence and one column per
#'   period, holding 1 for intervention, 0 for control and `NA` where clusters
#'   in that sequence are not part of the trial in that period (for example,
#'   clusters recruited at stage 2). Designs from [parallel_crt()] have a default
#'   in which sequences 1 and 2 are the original control and intervention
#'   clusters and sequences 3 and 4 the control and intervention clusters added
#'   at stage 2. If stage 2 may run for fewer periods than the schedule allows,
#'   give the longest schedule. Periods without data are treated as unobserved.
#' @param stage1 Integer vector of the consecutive periods forming stage 1.
#'   Defaults to period 1 for designs from [parallel_crt()].
#' @param strata Name of the column of the allocation data holding the
#'   randomisation stratum (for example `"country"`), or `NULL` if randomisation
#'   was unstratified. Give several names if randomisation was stratified on the
#'   cross-classification of several factors.
#' @param period_effects How to adjust for calendar time: `"fixed"` adds an
#'   additive effect for each period, `"none"` omits them, and `"auto"` (the
#'   default) uses `"fixed"` unless the proportion of sequences on intervention
#'   is the same in every period. See Details.
#' @param family One of `"binomial"`, `"poisson"` or `"gaussian"`.
#' @param w1 Stage 1 combination weight. The stage 2 weight is
#'   `sqrt(1 - w1^2)`.
#' @param c1 Stage 1 efficacy boundary for `|z1|`. Defaults to
#'   `qnorm(1 - alpha / 2) / w1`, the rule used by [adaptive_analysis()].
#' @param c2 Final boundary for `|Z|`. Defaults to the value from
#'   `calibrate_c2()`.
#' @param cac_plan Cluster autocorrelation used when it cannot be estimated at
#'   the interim: stage 1 has a single period, or the between-cluster variance
#'   is estimated at its floor.
#' @param sigma2 For the Gaussian family only: a fixed residual variance. By
#'   default this is estimated from the stage 1 within-cluster-period variation.
#' @param alpha Two-sided significance level. It is used for the descriptive
#'   intervals and for any boundary not supplied by `planned`.
#' @param control Named list overriding numerical settings: `cac_bounds`
#'   (`c(0, 0.99)`), `icc_floor` (`1e-4`), `p_bounds` (`c(0.01, 0.99)`), `gh_nodes`
#'   (`32`), `span` (`6`), `grid_min` (`80`) and `grid_max` (`2000`). The binomial
#'   and Poisson posterior grid has at least `grid_min` points and is refined, up
#'   to `grid_max`, when cluster-periods are so informative that the posterior
#'   would be narrower than the grid spacing.
#'
#' @details
#' The repeated confidence interval needs the critical value that controls
#' coverage over the whole range of `z1`. When the design used binding futility,
#' `c2` is calibrated on the continuation region only, so the interval is based
#' on the non-binding value from `calibrate_c2()` instead. The test decision
#' still uses the design's `c2`. The two can then disagree slightly.
#'
#' Stratification enters in three places: the intervention indicator is centred
#' within stratum and period, each stratum has its own intercept, and variance
#' components are estimated within strata. Strata too small to leave any degrees
#' of freedom for the variance components are ignored for that one step, with a
#' warning.
#'
#' Period effects are needed when treatment is confounded with calendar time, as
#' in a stepped-wedge design. In a balanced design they can be omitted, at the
#' cost of a conservative test if there are real secular trends. Effects for
#' stage 1 periods are estimated from the interim data and locked. Effects for
#' stage 2 periods cannot be known at the interim, so they are estimated from
#' stage 2 data at the final analysis, with the other nuisance values held at
#' their interim values. The stage 2 statistic stays valid given stage 1 because
#' the adaptation depends only on stage 1 data. Every period effect is
#' re-estimated under each hypothesised treatment effect.
#'
#' @return An `analysis_plan` object.
#' @seealso [interim_statistic()], [final_statistic()], [effect_estimates()]
#' @examples
#' # IMPALA-style four-period crossover, without a design object
#' plan <- analysis_plan(
#'   family = "binomial",
#'   schedule = rbind(c(0, 1, 0, 1), c(1, 0, 1, 0)),
#'   stage1 = 1:2,
#'   strata = "country",
#'   w1 = 0.7,
#'   cac_plan = 0.8
#' )
#' plan
#' @export
analysis_plan <- function(planned = NULL, schedule = NULL, stage1 = NULL,
                          strata = NULL,
                          period_effects = c("auto", "none", "fixed"),
                          family = NULL, w1 = NULL, c1 = NULL,
                          c2 = NULL, cac_plan = NULL, sigma2 = NULL,
                          alpha = 0.05, control = list()) {

  futility <- NULL
  period_effects <- match.arg(period_effects)

  if (!is.null(planned)) {
    if (!inherits(planned, "adaptive_results")) {
      stop("planned must be an adaptive_results object from adaptive_analysis()")
    }
    if (planned$explored) {
      stop("planned must be a single-design result. Use adaptive_analysis() ",
           "with the single stage 1 design chosen for the trial.")
    }
    params <- planned$raw$results$params
    if (is.null(params)) {
      stop("planned has no decision-rule results, so there are no weights or ",
           "boundaries to analyse with (", planned$raw$message, ").")
    }
    design <- planned$design
    design_family <- design$fixed_params$family
    if (!is.null(family) && !is.null(design_family) &&
        !identical(family, design_family)) {
      stop(sprintf("family = '%s' differs from the design's family '%s'.",
                   family, design_family))
    }
    family <- family %||% design_family %||% "gaussian"
    w1 <- w1 %||% params$w1_ref
    c1 <- c1 %||% params$efficacy_boundary
    c2 <- c2 %||% params$c2
    cac_plan <- cac_plan %||% design$fixed_params$cac
    futility <- params$futility

    if (identical(design$spec$design_type, "parallel")) {
      if (is.null(schedule)) {
        s1 <- design$stage1_grid[1, , drop = FALSE]
        schedule <- rbind(c(0, 0), c(1, 1))
        if (max(design$stage2_grid_fn(s1)$k2) > 0) {
          schedule <- rbind(schedule, c(NA, 0), c(NA, 1))
        }
      }
      stage1 <- stage1 %||% 1L
    }
  }

  if (is.null(family)) stop("family must be supplied when planned is NULL")
  fam <- .estimation_family(family)

  if (is.null(schedule) || is.null(stage1)) {
    stop("schedule and stage1 must be supplied for this design type")
  }
  schedule <- as.matrix(schedule)
  if (!is.numeric(schedule) && !all(is.na(schedule))) {
    stop("schedule must be a numeric matrix")
  }
  storage.mode(schedule) <- "double"
  if (!all(schedule %in% c(0, 1, NA))) {
    stop("schedule may contain only 0, 1 and NA")
  }
  stage1 <- as.integer(stage1)
  if (!length(stage1) || any(stage1 < 1L | stage1 > ncol(schedule)) ||
      any(diff(stage1) != 1L)) {
    stop("stage1 must be consecutive periods within the schedule")
  }
  if (length(stage1) >= ncol(schedule)) {
    stop("schedule must include at least one period after stage 1")
  }
  if (period_effects == "auto") {
    on_intervention <- colMeans(schedule, na.rm = TRUE)
    period_effects <- if (diff(range(on_intervention)) < 1e-8) "none" else "fixed"
  }

  if (is.null(w1)) stop("w1 must be supplied when planned is NULL")
  if (w1 <= 0 || w1 >= 1) stop("w1 must be strictly between 0 and 1")
  c1 <- c1 %||% (stats::qnorm(1 - alpha / 2) / w1)
  c2 <- c2 %||% calibrate_c2(w1, c1, alpha = alpha)
  c2_ci <- if (identical(futility, "binding")) {
    calibrate_c2(w1, c1, alpha = alpha)
  } else {
    c2
  }

  if (is.null(cac_plan)) stop("cac_plan must be supplied when planned is NULL")
  if (!is.null(sigma2) && family != "gaussian") {
    stop("sigma2 applies only to the gaussian family")
  }

  structure(
    list(
      planned = planned,
      family = family,
      fam = fam,
      schedule = schedule,
      periods = seq_len(ncol(schedule)),
      stage1 = stage1,
      strata = strata,
      period_effects = period_effects,
      w1 = w1,
      w2 = sqrt(1 - w1^2),
      c1 = c1,
      c2 = c2,
      c2_ci = c2_ci,
      futility = futility,
      cac_plan = cac_plan,
      sigma2 = sigma2,
      alpha = alpha,
      control = .estimation_control(control)
    ),
    class = "analysis_plan"
  )
}

#' @export
print.analysis_plan <- function(x, ...) {
  cat("Adaptive CRT analysis plan\n")
  cat(sprintf("  Family: %s\n", x$family))
  cat(sprintf("  Periods: %d (stage 1: %s)\n", ncol(x$schedule),
              paste(x$stage1, collapse = ", ")))
  cat(sprintf("  Sequences: %d\n", nrow(x$schedule)))
  cat(sprintf("  Strata: %s\n",
              if (is.null(x$strata)) "none" else paste(x$strata, collapse = " x ")))
  cat(sprintf("  Period effects: %s\n", x$period_effects))
  cat(sprintf("  Weights: w1 = %.4f, w2 = %.4f\n", x$w1, x$w2))
  cat(sprintf("  Boundaries: |z1| > %.3f, |Z| > %.3f\n", x$c1, x$c2))
  if (!isTRUE(all.equal(x$c2, x$c2_ci))) {
    cat(sprintf("  Repeated CI critical value: %.3f (non-binding)\n", x$c2_ci))
  }
  cat(sprintf("  Linked design: %s\n", if (is.null(x$planned)) "no" else "yes"))
  invisible(x)
}

# =============================================================================
# Preparing the allocation and the data
# =============================================================================

# Intervention indicator x_it and its design-based centring within
# stratum x period: x_tilde = x - mean_j x_jt over all clusters randomised and
# active in that period, so x_tilde does not depend on which data are observed.
.prepare_design <- function(plan, alloc) {
  alloc <- as.data.frame(alloc)
  needed <- c("cluster", "sequence", plan$strata)
  missing_cols <- setdiff(needed, names(alloc))
  if (length(missing_cols)) {
    stop("alloc is missing column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (anyNA(alloc[needed])) stop("alloc contains missing values")
  if (anyDuplicated(alloc$cluster)) stop("alloc must have one row per cluster")
  if (!all(alloc$sequence %in% seq_len(nrow(plan$schedule)))) {
    stop("alloc$sequence must index a row of the schedule")
  }

  key <- as.character(alloc$cluster)
  ord <- order(key)
  alloc <- alloc[ord, , drop = FALSE]
  clusters <- key[ord]
  N <- length(clusters)
  Tn <- ncol(plan$schedule)
  sequence <- as.integer(alloc$sequence)
  strata <- if (is.null(plan$strata)) {
    rep("all", N)
  } else {
    do.call(paste, c(lapply(alloc[plan$strata], as.character), sep = ":"))
  }

  S <- plan$schedule[sequence, , drop = FALSE]
  dimnames(S) <- NULL
  active <- !is.na(S)
  X <- S
  X[!active] <- 0

  XT <- matrix(0, N, Tn)
  for (s in unique(strata)) {
    in_stratum <- strata == s
    for (t in seq_len(Tn)) {
      on <- in_stratum & active[, t]
      if (any(on)) XT[on, t] <- X[on, t] - mean(X[on, t])
    }
  }

  contrast <- tapply(rowSums(abs(XT)), strata, sum)
  if (any(contrast == 0) && length(contrast) > 1L) {
    warning("No treatment contrast within stratum/strata: ",
            paste(names(contrast)[contrast == 0], collapse = ", "),
            ". These clusters contribute nothing to the test.", call. = FALSE)
  }

  list(clusters = clusters, strata = strata, sequence = sequence, N = N,
       Tn = Tn, active = active, X = X, XT = XT, alloc = alloc[needed])
}

# Individual-level records to cluster-period sums. Records with a missing
# outcome are excluded.
.aggregate_individual <- function(data, family, cell) {
  keep <- !is.na(data$outcome)
  outcome <- data$outcome[keep]
  cell <- cell[keep]

  if (family == "binomial" && !all(outcome %in% c(0, 1))) {
    stop("outcome must be 0, 1 or NA for the binomial family")
  }
  if (family == "poisson" && any(outcome < 0)) {
    stop("outcome must be non-negative for the poisson family")
  }

  m <- rep(1, length(outcome))
  if (family == "poisson" && "exposure" %in% names(data)) {
    m <- data$exposure[keep]
    if (anyNA(m) || any(m <= 0)) {
      stop("exposure must be positive for every record with an outcome")
    }
  }

  centre <- if (length(outcome)) mean(outcome) else 0
  shifted <- outcome - centre
  sums <- rowsum(cbind(y = outcome, m = m, n = 1, shift = shifted,
                       shift2 = shifted^2), cell)
  list(
    cell = as.integer(rownames(sums)),
    y = sums[, "y"],
    m = sums[, "m"],
    ss = sums[, "shift2"] - sums[, "shift"]^2 / sums[, "n"],
    df = sums[, "n"] - 1
  )
}

# Cluster-period records, with validation.
.read_cluster_period <- function(data, family, cell) {
  if (anyDuplicated(cell)) {
    stop("data must have one row per cluster-period")
  }
  if (family == "gaussian" && "ybar" %in% names(data) &&
      !"y" %in% names(data)) {
    data$y <- data$ybar * data$m
  }
  needed <- c("y", "m")
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("data needs an 'outcome' column or column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  m <- data$m
  m[is.na(m)] <- 0
  y <- data$y
  keep <- m > 0
  if (any(m < 0)) stop("m must be non-negative")
  if (anyNA(y[keep])) stop("y is missing for a cluster-period with m > 0")
  if (family == "binomial" && any(y[keep] < 0 | y[keep] > m[keep])) {
    stop("y must lie between 0 and m for the binomial family")
  }
  if (family == "poisson" && any(y[keep] < 0)) {
    stop("y must be non-negative for the poisson family")
  }

  ss <- df <- NULL
  if (family == "gaussian" && "s2" %in% names(data)) {
    df <- pmax(m[keep] - 1, 0)
    ss <- data$s2[keep] * df
    ss[df == 0] <- 0
    if (anyNA(ss)) stop("s2 is missing for a cluster-period with m > 1")
  }
  list(cell = cell[keep], y = y[keep], m = m[keep], ss = ss, df = df)
}

# Cluster x period matrices of summed outcomes (Y) and denominators (M).
.prepare_data <- function(plan, design, data) {
  data <- as.data.frame(data)
  missing_cols <- setdiff(c("cluster", "period"), names(data))
  if (length(missing_cols)) {
    stop("data is missing column(s): ", paste(missing_cols, collapse = ", "))
  }
  i <- match(as.character(data$cluster), design$clusters)
  if (anyNA(i)) {
    stop("data contains cluster(s) not in the allocation: ",
         paste(utils::head(unique(data$cluster[is.na(i)]), 5), collapse = ", "))
  }
  if (!all(data$period %in% plan$periods)) {
    stop("data$period must be an integer between 1 and ", ncol(plan$schedule))
  }
  t <- as.integer(data$period)
  N <- design$N
  Tn <- design$Tn
  cell <- (t - 1L) * N + i

  cells <- if ("outcome" %in% names(data)) {
    .aggregate_individual(data, plan$family, cell)
  } else {
    .read_cluster_period(data, plan$family, cell)
  }

  Y <- M <- matrix(0, N, Tn)
  Y[cells$cell] <- cells$y
  M[cells$cell] <- cells$m
  SS <- DF <- NULL
  if (!is.null(cells$ss)) {
    SS <- DF <- matrix(0, N, Tn)
    SS[cells$cell] <- cells$ss
    DF[cells$cell] <- cells$df
  }
  obs <- M > 0
  if (any(obs & !design$active)) {
    stop("data contains outcomes for a cluster-period the schedule marks as ",
         "not part of the trial (NA)")
  }

  period_cells <- NULL
  if (plan$period_effects == "fixed" && any(obs)) {
    keep <- which(obs)
    period_cells <- stats::aggregate(
      cbind(m = M[keep], y = Y[keep]),
      list(stratum = design$strata[row(M)[keep]], period = col(M)[keep],
           x = design$X[keep]),
      sum
    )
  }
  list(Y = Y, M = M, SS = SS, DF = DF, obs = obs, cells = period_cells)
}

# Warn if the observed stage 1 trial differs from the planned stage 1 design.
.check_planned_clusters <- function(plan, design) {
  if (is.null(plan$planned)) return(invisible(NULL))
  planned_design <- plan$planned$design
  formula <- planned_design$spec$resources$clusters_s1
  if (is.null(formula)) return(invisible(NULL))
  env <- c(as.list(planned_design$stage1_grid[1, , drop = FALSE]),
           list(n_arms = planned_design$spec$n_arms))
  expected <- tryCatch(eval_resource(formula, env), error = function(e) NA_real_)
  actual <- sum(rowSums(design$active[, plan$stage1, drop = FALSE]) > 0)
  if (!is.na(expected) && expected != actual) {
    warning(sprintf(
      "The allocation has %d stage 1 clusters but the design planned %g.",
      actual, expected
    ), call. = FALSE)
  }
  invisible(NULL)
}

# =============================================================================
# Interim nuisance estimation (stage 1 data only)
# =============================================================================

.estimate_nuisance <- function(plan, design, d) {
  fam <- plan$fam
  control <- plan$control
  s1 <- plan$stage1

  Y <- d$Y[, s1, drop = FALSE]
  M <- d$M[, s1, drop = FALSE]
  X <- design$X[, s1, drop = FALSE]
  obs <- d$obs[, s1, drop = FALSE]

  # (a) Control-arm mean. Used only to choose the stage 2 design.
  control_cells <- obs & X == 0
  m_control <- sum(M[control_cells])
  if (m_control <= 0) stop("No stage 1 outcomes observed under control.")
  baseline <- fam$clamp_mean(sum(Y[control_cells]) / m_control, m_control,
                             control)

  sigma2 <- NA_real_
  if (plan$family == "gaussian") {
    sigma2 <- plan$sigma2
    if (is.null(sigma2)) {
      if (is.null(d$SS)) {
        stop("The gaussian family needs individual-level data, an 's2' column ",
             "in the cluster-period data, or sigma2 in analysis_plan().")
      }
      df_within <- sum(d$DF[, s1])
      if (df_within < 1) {
        stop("No within-cluster-period degrees of freedom to estimate sigma2.")
      }
      sigma2 <- sum(d$SS[, s1]) / df_within
    }
  }

  # (b) Working responses and their sampling variances, for clusters observed
  # in every stage 1 period, centred within cells of stratum x stage 1 treatment
  # pattern x period. Sequences that do not differ during stage 1 share a cell,
  # which matters for designs with one cluster per sequence.
  complete <- rowSums(obs) == length(s1)
  n_complete <- sum(complete)
  pattern <- apply(X, 1, paste, collapse = "")[complete]
  group <- paste(design$strata[complete], pattern, sep = "\r")
  df <- n_complete - length(unique(group))
  if (df < 1 && length(unique(design$strata)) > 1L) {
    group <- pattern
    df <- n_complete - length(unique(group))
    if (df >= 1) {
      warning("Strata are too small to estimate variance components within ",
              "them. Ignoring strata for this step only.", call. = FALSE)
    }
  }
  if (df < 1) stop("No degrees of freedom for variance estimation.")

  wk <- fam$working(Y[complete, , drop = FALSE], M[complete, , drop = FALSE],
                    sigma2)
  ell <- matrix(wk$ell, n_complete, length(s1))
  w <- matrix(wk$w, n_complete, length(s1))
  n_g <- stats::ave(rep(1, n_complete), group, FUN = length)
  r <- apply(ell, 2, function(v) v - stats::ave(v, group))

  A_t <- (colSums(r^2) - colSums((n_g - 1) / n_g * w)) / df
  A_raw <- mean(A_t)
  s_raw <- if (length(s1) >= 2L) {
    mean(vapply(seq_len(length(s1) - 1L), function(t) {
      sum(r[, t] * r[, t + 1L]) / df
    }, numeric(1)))
  } else {
    NA_real_
  }

  # (c) Truncation: ICC >= icc_floor, CAC within cac_bounds.
  scale <- fam$icc_scale(baseline, sigma2)
  A_min <- control$icc_floor / (1 - control$icc_floor) / scale
  A <- max(A_raw, A_min)
  cac <- if (A_raw <= A_min || is.na(s_raw)) {
    plan$cac_plan
  } else {
    .clamp(s_raw / A, control$cac_bounds)
  }
  q <- A * scale
  icc <- q / (1 + q)

  # (d) Stratum summaries used to profile the intercepts under any beta0:
  # pooled stage 1 mean (both arms) and share of the denominator on intervention.
  m_all <- ifelse(obs, M, 0)
  strata_ids <- factor(design$strata)
  m_s <- tapply(rowSums(m_all), strata_ids, sum)
  y_s <- tapply(rowSums(ifelse(obs, Y, 0)), strata_ids, sum)
  x_s <- tapply(rowSums(m_all * X), strata_ids, sum)
  seen <- which(m_s > 0)
  strata_tab <- data.frame(
    stratum = names(m_s)[seen],
    mean = vapply(seen, function(j) {
      fam$clamp_mean(y_s[[j]] / m_s[[j]], m_s[[j]], control)
    }, numeric(1)),
    share1 = as.numeric(x_s[seen] / m_s[seen]),
    stringsAsFactors = FALSE
  )

  gh <- NULL
  if (plan$family == "binomial") {
    g <- statmod::gauss.quad(control$gh_nodes, kind = "hermite")
    gh <- list(nodes = g$nodes * sqrt(2), weights = g$weights / sqrt(pi))
  }

  nuis <- list(
    A = A, cac = cac, icc = icc, baseline = baseline, sigma2 = sigma2,
    strata = strata_tab, gh = gh, cells = d$cells,
    grid = if (fam$grid) .posterior_grid(A, cac, control) else NULL,
    raw = list(A = A_raw,
               cac = if (A_raw > 0) s_raw / A_raw else NA_real_,
               df = df, n_clusters = n_complete)
  )
  nuis$eta <- .intercepts(plan, nuis, 0)
  nuis
}

# Stratum intercepts on the cluster-specific scale under hypothesised beta,
# solving share0 * E[mu(eta + b)] + share1 * E[mu(eta + beta + b)] = pooled mean.
.intercepts <- function(plan, nuis, beta) {
  s <- nuis$strata
  stats::setNames(
    vapply(seq_len(nrow(s)), function(j) {
      plan$fam$intercept(s$mean[j], s$share1[j], beta, nuis$A, nuis$gh)
    }, numeric(1)),
    s$stratum
  )
}

# Newton solve of X' (m * E[mu(X theta + offset + b)]) = target.
.marginal_fit <- function(fam, X, offset, target, m, A, gh, start,
                          maxit = 50L, tol = 1e-10) {
  theta <- start
  size <- pmax(drop(crossprod(X, m)), 1)
  for (it in seq_len(maxit)) {
    eta <- drop(X %*% theta) + offset
    score <- target - drop(crossprod(X, m * fam$marginal(eta, A, gh)))
    if (max(abs(score) / size) < tol) break
    W <- m * fam$dmarginal(eta, A, gh)
    H <- crossprod(X, X * W)
    step <- solve(H + diag(1e-10, ncol(X)), score)
    if (fam$grid && max(abs(step)) > 2) step <- step * 2 / max(abs(step))
    theta <- theta + step
  }
  theta
}

# Stratum intercepts and period effects under hypothesised beta. Stage 1 effects
# use the data locked at the interim. Stage 2 effects use the data in `d`, with
# the stratum intercepts held at their interim values.
.linear_pars <- function(plan, nuis, beta, d, Tn) {
  eta <- .intercepts(plan, nuis, beta)
  gamma <- numeric(Tn)
  if (plan$period_effects == "none") return(list(eta = eta, gamma = gamma))

  fam <- plan$fam
  control <- plan$control
  target_mean <- function(y, m) m * fam$clamp_mean(y / m, m, control)
  later <- plan$stage1[-1]

  if (length(later)) {
    cells <- nuis$cells
    cells <- cells[cells$period %in% plan$stage1, , drop = FALSE]
    s <- match(cells$stratum, nuis$strata$stratum)
    X <- cbind(diag(nrow(nuis$strata))[s, , drop = FALSE],
               outer(cells$period, later, "==") + 0)
    m_s <- tapply(cells$m, factor(s, seq_len(nrow(nuis$strata))), sum)
    target <- c(as.numeric(m_s) * nuis$strata$mean,
                vapply(later, function(t) {
                  on <- cells$period == t
                  if (!any(on)) return(0)
                  target_mean(sum(cells$y[on]), sum(cells$m[on]))
                }, numeric(1)))
    theta <- .marginal_fit(fam, X, beta * cells$x, target, cells$m, nuis$A,
                           nuis$gh, c(eta, numeric(length(later))))
    eta <- stats::setNames(theta[seq_len(nrow(nuis$strata))], names(eta))
    gamma[later] <- theta[-seq_len(nrow(nuis$strata))]
  }

  last <- gamma[max(plan$stage1)]
  if (!is.null(d$cells)) {
    for (t in setdiff(seq_len(Tn), plan$stage1)) {
      cells <- d$cells[d$cells$period == t, , drop = FALSE]
      s <- match(cells$stratum, names(eta))
      cells <- cells[!is.na(s), , drop = FALSE]
      s <- s[!is.na(s)]
      if (!nrow(cells)) next
      gamma[t] <- .marginal_fit(
        fam, matrix(1, nrow(cells), 1), unname(eta[s]) + beta * cells$x,
        target_mean(sum(cells$y), sum(cells$m)), cells$m, nuis$A, nuis$gh, last
      )
      last <- gamma[t]
    }
  }
  list(eta = eta, gamma = gamma)
}

# =============================================================================
# Posterior mean of the cluster-period mean under the AR(1) random effect
# =============================================================================

.posterior_grid <- function(A, cac, control) {
  span <- control$span
  ngrid <- max(as.integer(control$grid_min),
               ceiling(3 * 2 * span / sqrt(1 - cac^2)) + 1L)
  b <- seq(-span, span, length.out = ngrid) * sqrt(A)
  prior <- stats::dnorm(b, 0, sqrt(A))
  prior <- prior / sum(prior)
  K <- outer(b, cac * b, function(u, v) {
    stats::dnorm(u, v, sqrt(A * (1 - cac^2)))
  })
  # K[j, k] = P(b_t = b_j | b_{t-1} = b_k)
  K <- sweep(K, 2, colSums(K), "/")
  list(b = b, prior = prior, K = K)
}

# Grid spacing must not exceed the narrowest posterior standard deviation in
# the data, or a very informative cluster-period falls between grid points.
.resolve_grid <- function(plan, nuis, d) {
  grid <- nuis$grid
  if (!any(d$obs)) return(grid)
  control <- plan$control
  info <- max(plan$fam$info(d$Y[d$obs], d$M[d$obs]))
  precision <- info + 1 / (nuis$A * (1 - nuis$cac^2))
  needed <- ceiling(2 * control$span * sqrt(nuis$A * precision)) + 1L
  needed <- min(needed, as.integer(control$grid_max))
  if (needed <= length(grid$b)) return(grid)
  control$grid_min <- needed
  .posterior_grid(nuis$A, nuis$cac, control)
}

# Forward-backward on the grid, for families without a closed form.
.posterior_mean_grid <- function(fam, Y, M, ETA, obs, grid) {
  ncl <- nrow(Y)
  Tn <- ncol(Y)
  b <- grid$b
  K <- grid$K
  ngrid <- length(b)

  G <- lapply(seq_len(Tn), function(t) {
    L <- matrix(1, ncl, ngrid)              # unobserved: uninformative
    o <- obs[, t]
    if (any(o)) {
      H <- outer(ETA[o, t], b, "+")
      LL <- fam$loglik(H, Y[o, t], M[o, t])
      L[o, ] <- exp(LL - apply(LL, 1, max))
    }
    L
  })
  nrm <- function(x) x / pmax(rowSums(x), .Machine$double.xmin)

  fwd <- bwd <- vector("list", Tn)
  fwd[[1]] <- nrm(sweep(G[[1]], 2, grid$prior, "*"))
  for (t in seq_len(Tn)[-1]) fwd[[t]] <- nrm(G[[t]] * (fwd[[t - 1]] %*% t(K)))
  bwd[[Tn]] <- matrix(1, ncl, ngrid)
  for (t in rev(seq_len(Tn - 1))) {
    bwd[[t]] <- nrm((G[[t + 1]] * bwd[[t + 1]]) %*% K)
  }

  out <- vapply(seq_len(Tn), function(t) {
    g <- fwd[[t]] * bwd[[t]]
    rowSums(g * fam$linkinv(outer(ETA[, t], b, "+"))) /
      pmax(rowSums(g), .Machine$double.xmin)
  }, numeric(ncl))
  matrix(out, ncl, Tn)
}

# Exact conditional mean for the Gaussian family.
.posterior_mean_gaussian <- function(Y, M, ETA, obs, A, cac, sigma2) {
  Sigma <- A * cac^abs(outer(seq_len(ncol(Y)), seq_len(ncol(Y)), "-"))
  out <- ETA
  for (i in seq_len(nrow(Y))) {
    o <- which(obs[i, ])
    if (!length(o)) next
    resid <- Y[i, o] / M[i, o] - ETA[i, o]
    S <- Sigma[o, o, drop = FALSE]
    V <- S + diag(sigma2 / M[i, o], length(o))
    out[i, o] <- ETA[i, o] + drop(S %*% solve(V, resid))
  }
  out
}

# =============================================================================
# Score contributions and stage-wise statistics
# =============================================================================

# Cluster score contributions U_i(beta), from the data in `d`.
.score_contributions <- function(plan, design, d, nuis, beta = 0) {
  pars <- .linear_pars(plan, nuis, beta, d, design$Tn)
  eta <- pars$eta
  has_data <- rowSums(d$obs) > 0
  unknown <- setdiff(unique(design$strata[has_data]), names(eta))
  if (length(unknown)) {
    stop("No stage 1 data for stratum/strata: ",
         paste(unknown, collapse = ", "),
         ". Its intercept cannot be estimated at the interim.")
  }
  stratum_eta <- unname(eta[design$strata])
  stratum_eta[is.na(stratum_eta)] <- 0
  ETA <- matrix(stratum_eta, design$N, design$Tn) +
    matrix(pars$gamma, design$N, design$Tn, byrow = TRUE) + beta * design$X

  mu <- if (plan$fam$grid) {
    .posterior_mean_grid(plan$fam, d$Y, d$M, ETA, d$obs,
                         .resolve_grid(plan, nuis, d))
  } else {
    .posterior_mean_gaussian(d$Y, d$M, ETA, d$obs, nuis$A, nuis$cac,
                             nuis$sigma2)
  }
  stats::setNames(rowSums(design$XT * (d$Y - d$M * mu) * d$obs),
                  design$clusters)
}

.self_normalised <- function(u) {
  v <- sum(u^2)
  if (v <= 0) return(NA_real_)
  sum(u) / sqrt(v)
}

# Place interim contributions on the final cluster list. Clusters recruited
# after the interim contribute zero.
.align_interim_score <- function(U1, clusters) {
  out <- stats::setNames(numeric(length(clusters)), clusters)
  out[names(U1)] <- U1
  out
}

# Allocation at final lock: the interim allocation, plus any clusters recruited
# at stage 2. Interim clusters and their stage 1 centring cannot change.
.final_design <- function(interim, alloc) {
  if (is.null(alloc)) return(interim$design)
  plan <- interim$plan
  old <- interim$design
  new <- .prepare_design(plan, alloc)
  pos <- match(old$clusters, new$clusters)
  if (anyNA(pos)) {
    stop("alloc must contain every cluster in the interim allocation")
  }
  added <- setdiff(seq_len(new$N), pos)
  if (length(added) &&
      any(new$active[added, plan$stage1, drop = FALSE])) {
    stop("Clusters added after the interim must be NA in the schedule for all ",
         "stage 1 periods")
  }
  s1 <- plan$stage1
  if (!identical(new$strata[pos], old$strata) ||
      !identical(new$active[pos, s1, drop = FALSE],
                 old$active[, s1, drop = FALSE]) ||
      !isTRUE(all.equal(new$XT[pos, s1, drop = FALSE],
                        old$XT[, s1, drop = FALSE]))) {
    stop("alloc or schedule changes the stratum or stage 1 treatment of a ",
         "cluster used at the interim")
  }
  new
}

# Statistics at beta. `fin` is list(design, data) or NULL at the interim.
.stage_stats <- function(interim, fin, beta, which = c("z1", "z2", "Z", "z_all")) {
  plan <- interim$plan
  out <- list()
  needs_u1 <- any(c("z1", "z2", "Z") %in% which)
  needs_uf <- any(c("z2", "Z", "z_all") %in% which)

  if (needs_uf && is.null(fin)) {
    stop("Final data are required for z2, Z and z_all")
  }
  if (needs_u1) {
    U1 <- .score_contributions(plan, interim$design, interim$data, interim$nuis,
                               beta)
    out$z1 <- .self_normalised(U1)
  }
  if (needs_uf) {
    Uf <- .score_contributions(plan, fin$design, fin$data, interim$nuis, beta)
    out$z_all <- .self_normalised(Uf)
    if (needs_u1) {
      U2 <- Uf - .align_interim_score(U1, fin$design$clusters)
      out$z2 <- .self_normalised(U2)
      out$Z <- plan$w1 * out$z1 + plan$w2 * out$z2
    }
  }
  out
}

#' Interim analysis statistic for an adaptive cluster randomised trial
#'
#' Estimates the nuisance parameters from stage 1 data and computes the stage 1
#' score statistic. The nuisance values are locked here and reused at the final
#' analysis.
#'
#' @param plan An [analysis_plan()] object.
#' @param data Data locked for the interim, in one of two forms. Individual-level
#'   data has columns `cluster`, `period` and `outcome` (`NA` for a missing
#'   outcome), and for Poisson optionally `exposure`. Cluster-period data has
#'   columns `cluster`, `period`, `y` and `m`: `y` is the event count or the sum
#'   of outcomes and `m` is the number of participants with an observed outcome
#'   or the exposure. For Gaussian cluster-period data, `ybar` can replace `y`,
#'   and a within-cluster-period variance column `s2` is needed unless `sigma2`
#'   is set in the plan. Rows for periods after stage 1 are ignored.
#' @param alloc Data frame with one row per randomised cluster: `cluster`,
#'   `sequence` (a row of the schedule) and the strata column, if any.
#'
#' @details
#' Exposure in Poisson data must use the same units as the design's `baseline`,
#' because the ICC is defined on that rate scale. For the Gaussian family
#' `theta_hat$icc` is `A / (A + sigma2)`.
#'
#' @return An `interim_statistic` object, with elements `z1`, `U1` (cluster score
#'   contributions), `nuis` (locked nuisance values) and `theta_hat` (`icc`, `cac`
#'   and `baseline`, ready for [interim_analysis()]).
#' @seealso [interim_decision()], [final_statistic()], [effect_estimates()]
#' @examples
#' plan <- analysis_plan(family = "binomial",
#'                       schedule = rbind(c(0, 1, 0, 1), c(1, 0, 1, 0)),
#'                       stage1 = 1:2, w1 = 0.7, cac_plan = 0.8)
#' alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
#' set.seed(1)
#' effect <- rnorm(20, 0, 0.3)
#' data_all <- expand.grid(cluster = 1:20, period = 1:4)
#' data_all$m <- 40
#' data_all$y <- rbinom(80, 40, plogis(-1 + effect[data_all$cluster]))
#' interim <- interim_statistic(plan, data_all[data_all$period <= 2, ], alloc)
#' interim
#' interim$theta_hat
#' @export
interim_statistic <- function(plan, data, alloc) {
  if (!inherits(plan, "analysis_plan")) {
    stop("plan must be an analysis_plan object")
  }
  design <- .prepare_design(plan, alloc)
  .check_planned_clusters(plan, design)

  data <- as.data.frame(data)
  if ("period" %in% names(data)) {
    in_stage1 <- data$period %in% plan$stage1
    if (!all(in_stage1)) {
      message(sprintf("Ignoring %d row(s) from periods after stage 1.",
                      sum(!in_stage1)))
      data <- data[in_stage1, , drop = FALSE]
    }
  }
  d <- .prepare_data(plan, design, data)
  nuis <- .estimate_nuisance(plan, design, d)
  U1 <- .score_contributions(plan, design, d, nuis, 0)

  structure(
    list(
      z1 = .self_normalised(U1),
      U1 = U1,
      nuis = nuis,
      theta_hat = nuis[c("icc", "cac", "baseline")],
      plan = plan,
      design = design,
      data = d
    ),
    class = "interim_statistic"
  )
}

#' @export
print.interim_statistic <- function(x, ...) {
  cat("Interim statistic\n")
  cat(sprintf("  z1 = %.3f (efficacy boundary: |z1| > %.3f)\n", x$z1, x$plan$c1))
  cat(sprintf("  ICC = %.4f, CAC = %.3f, baseline = %.4g\n",
              x$nuis$icc, x$nuis$cac, x$nuis$baseline))
  cat(sprintf("  Clusters in variance estimation: %d (df = %d)\n",
              x$nuis$raw$n_clusters, x$nuis$raw$df))
  invisible(x)
}

#' Decide what to do at the interim analysis
#'
#' Passes an [interim_statistic()] to [interim_analysis()], using the design the
#' plan was built from.
#'
#' @param interim An `interim_statistic` object whose plan was built with
#'   `planned`.
#' @param ... Further arguments for [interim_analysis()], for example `delta`,
#'   `recalibrate` or `verbose`. `target_power` defaults to the planned value.
#' @return An `interim_result` object.
#' @examples
#' \donttest{
#' design <- parallel_crt(icc = 0.05, cac = 0.8, delta = 0.25, k1 = 14, m1 = 30,
#'                        k2 = 0:2, m2 = c(20, 30), rho = 30)
#' planned <- adaptive_analysis(design, target_power = 0.8, tol = 0.02,
#'                              method = "lambda", verbose = FALSE)
#' plan <- analysis_plan(planned)
#' alloc <- data.frame(cluster = 1:28, sequence = rep(1:2, each = 14))
#' stage1 <- data.frame(cluster = 1:28, period = 1, m = 30,
#'                      y = rnorm(28, rep(c(0, 0.25), each = 14), 0.2) * 30, s2 = 1)
#' interim <- interim_statistic(plan, stage1, alloc)
#' interim_decision(interim, verbose = FALSE)
#' }
#' @export
interim_decision <- function(interim, ...) {
  if (!inherits(interim, "interim_statistic")) {
    stop("interim must be an interim_statistic object")
  }
  planned <- interim$plan$planned
  if (is.null(planned)) {
    stop("This analysis plan is not linked to a design. Build it with ",
         "analysis_plan(planned = ...) or call interim_analysis() directly.")
  }
  args <- utils::modifyList(list(target_power = planned$target_power), list(...))
  do.call(interim_analysis,
          c(list(planned = planned, z1_obs = interim$z1,
                 theta_hat = interim$theta_hat), args))
}

#' Final combination test for an adaptive cluster randomised trial
#'
#' Computes the stage 2 score increment and the combination statistic, reusing
#' `z1` and the nuisance values locked at the interim.
#'
#' @param interim An [interim_statistic()] object.
#' @param data All data at final lock, including any stage 1 outcomes
#'   ascertained after the interim lock. Same format as in [interim_statistic()].
#' @param alloc Allocation at final lock. Needed only if clusters were added at
#'   stage 2. The clusters used at the interim must keep their stage 1 treatment
#'   and stratum, and added clusters must use a sequence that is `NA` for all
#'   stage 1 periods.
#' @param schedule Schedule at final lock. Needed only if the interim decision
#'   changed the stage 2 schedule, for example the number of stage 2 periods or
#'   the rollout of a stepped-wedge design. Stage 1 treatment must be unchanged.
#'
#' @return A `final_statistic` object with elements `z1`, `z2`, `Z`, `reject`
#'   (`|Z| > c2`) and `U2` (cluster score increments).
#' @seealso [effect_estimates()], [stage_statistics()]
#' @examples
#' plan <- analysis_plan(family = "binomial",
#'                       schedule = rbind(c(0, 1, 0, 1), c(1, 0, 1, 0)),
#'                       stage1 = 1:2, w1 = 0.7, cac_plan = 0.8)
#' alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
#' set.seed(1)
#' effect <- rnorm(20, 0, 0.3)
#' data_all <- expand.grid(cluster = 1:20, period = 1:4)
#' data_all$m <- 40
#' data_all$y <- rbinom(80, 40, plogis(-1 + effect[data_all$cluster]))
#' interim <- interim_statistic(plan, data_all[data_all$period <= 2, ], alloc)
#' final_statistic(interim, data_all)
#' @export
final_statistic <- function(interim, data, alloc = NULL, schedule = NULL) {
  if (!inherits(interim, "interim_statistic")) {
    stop("interim must be an interim_statistic object")
  }
  plan <- interim$plan
  if (!is.null(schedule)) {
    schedule <- as.matrix(schedule)
    storage.mode(schedule) <- "double"
    if (!all(schedule %in% c(0, 1, NA)) ||
        ncol(schedule) <= max(plan$stage1)) {
      stop("schedule may contain only 0, 1 and NA, and must include a period ",
           "after stage 1")
    }
    plan$schedule <- schedule
    plan$periods <- seq_len(ncol(schedule))
    interim$plan <- plan
    if (is.null(alloc)) alloc <- interim$design$alloc
  }
  design <- .final_design(interim, alloc)
  d <- .prepare_data(plan, design, data)
  Uf <- .score_contributions(plan, design, d, interim$nuis, 0)
  U2 <- Uf - .align_interim_score(interim$U1, design$clusters)
  z2 <- .self_normalised(U2)
  Z <- plan$w1 * interim$z1 + plan$w2 * z2

  structure(
    list(
      z1 = interim$z1,
      z2 = z2,
      Z = Z,
      c2 = plan$c2,
      reject = abs(Z) > plan$c2,
      U2 = U2,
      interim = interim,
      design = design,
      data = d
    ),
    class = "final_statistic"
  )
}

#' @export
print.final_statistic <- function(x, ...) {
  cat("Final combination test\n")
  cat(sprintf("  z1 = %.3f, z2 = %.3f\n", x$z1, x$z2))
  cat(sprintf("  Z = %.3f (boundary: |Z| > %.3f)\n", x$Z, x$c2))
  cat(sprintf("  Reject null hypothesis: %s\n", if (isTRUE(x$reject)) "yes" else "no"))
  invisible(x)
}

#' Stage-wise statistics at a hypothesised treatment effect
#'
#' Evaluates the stage-wise statistics at `beta` on the link scale, using the
#' interim nuisance values with the stratum intercepts re-profiled under `beta`.
#' At `beta = 0` the values equal those from [interim_statistic()] and
#' [final_statistic()]. Each statistic is expected to decrease in `beta`.
#'
#' @param object An [interim_statistic()] or [final_statistic()] object.
#' @param beta Hypothesised treatment effect on the link scale.
#' @return A list with `z1` and, for a `final_statistic`, `z2`, `Z` and `z_all`
#'   (the statistic from all data, ignoring the adaptation).
#' @examples
#' plan <- analysis_plan(family = "binomial",
#'                       schedule = rbind(c(0, 1, 0, 1), c(1, 0, 1, 0)),
#'                       stage1 = 1:2, w1 = 0.7, cac_plan = 0.8)
#' alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
#' set.seed(1)
#' effect <- rnorm(20, 0, 0.3)
#' data_all <- expand.grid(cluster = 1:20, period = 1:4)
#' data_all$m <- 40
#' data_all$y <- rbinom(80, 40, plogis(-1 + effect[data_all$cluster]))
#' interim <- interim_statistic(plan, data_all[data_all$period <= 2, ], alloc)
#' final <- final_statistic(interim, data_all)
#' stage_statistics(final, beta = 0.2)
#' @export
stage_statistics <- function(object, beta = 0) {
  if (inherits(object, "final_statistic")) {
    .stage_stats(object$interim, object[c("design", "data")], beta)
  } else if (inherits(object, "interim_statistic")) {
    .stage_stats(object, NULL, beta, "z1")
  } else {
    stop("object must be an interim_statistic or final_statistic object")
  }
}

# =============================================================================
# Effect estimation
# =============================================================================

# Root of a decreasing function, expanding the search from `interval` up to
# +/- limit. Returns NA if no root can be bracketed.
.decreasing_root <- function(g, interval, limit, tol) {
  value <- function(b) tryCatch(g(b), error = function(e) NA_real_)
  lo <- interval[1]
  hi <- interval[2]
  f_lo <- value(lo)
  f_hi <- value(hi)
  width <- hi - lo
  while (is.finite(f_lo) && f_lo <= 0 && lo > -limit) {
    width <- width * 2
    lo <- max(lo - width, -limit)
    f_lo <- value(lo)
  }
  width <- hi - lo
  while (is.finite(f_hi) && f_hi >= 0 && hi < limit) {
    width <- width * 2
    hi <- min(hi + width, limit)
    f_hi <- value(hi)
  }
  if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo <= 0 || f_hi >= 0) {
    return(NA_real_)
  }
  stats::uniroot(g, c(lo, hi), f.lower = f_lo, f.upper = f_hi, tol = tol)$root
}

# Estimate, interval and standard error from a decreasing statistic f(beta0):
#   estimate: f(estimate) = 0
#   interval: {beta0 : |f(beta0)| <= crit}
#   se:       -1 / f'(estimate), so that f(beta0) ~ (estimate - beta0) / se
.invert_statistic <- function(f, crit, scale) {
  interval <- c(-3, 3) * scale
  limit <- 30 * scale
  tol <- 1e-7 * scale
  h <- 1e-3 * scale
  root <- function(target) {
    .decreasing_root(function(b) f(b) - target, interval, limit, tol)
  }
  est <- root(0)
  lower <- if (is.finite(crit)) root(crit) else -Inf
  upper <- if (is.finite(crit)) root(-crit) else Inf
  if (is.na(lower)) lower <- -Inf
  if (is.na(upper)) upper <- Inf
  se <- if (is.finite(est)) 2 * h / (f(est - h) - f(est + h)) else NA_real_
  c(estimate = est, lower = lower, upper = upper, se = se)
}

#' Treatment effect estimates and repeated confidence intervals
#'
#' Inverts the stage-wise statistics to give point estimates and intervals for
#' the treatment effect.
#'
#' @param object A [final_statistic()] object, or an [interim_statistic()] object
#'   if the trial stopped at the interim.
#' @param check If `TRUE`, warn when the primary statistic is not monotone in the
#'   hypothesised effect across and beyond its interval.
#'
#' @details
#' The rows of the result are:
#' \describe{
#'   \item{`adaptive`}{Primary: the repeated confidence interval at the stage the
#'   trial stopped, using `|Z| <= c2` at the final analysis or `|z1| <= c1` at the
#'   interim.}
#'   \item{`stage1`, `stage2`}{Stage-wise estimates with unadjusted intervals
#'   (descriptive). Final analysis only.}
#'   \item{`unadjusted`}{All data analysed as a fixed design, ignoring the
#'   adaptation.}
#' }
#'
#' `estimate`, `lower`, `upper` and `se` are on the link scale: the log odds
#' ratio (binomial), log rate ratio (Poisson) or mean difference (Gaussian). The
#' `ratio` columns are the exponentiated values and are `NA` for the Gaussian
#' family.
#'
#' The `difference` columns are the marginal risk or rate difference implied by
#' the link-scale effect, with strata weighted by their share of participants
#' (or exposure) with an observed outcome. Interval limits are transformed
#' endpoints and the standard error uses the delta method.
#'
#' Self-normalised statistics are bounded, so with few clusters an interval
#' limit may not exist. It is then returned as infinite with a warning.
#'
#' @return A data frame of class `effect_estimates`, one row per analysis.
#' @examples
#' plan <- analysis_plan(family = "binomial",
#'                       schedule = rbind(c(0, 1, 0, 1), c(1, 0, 1, 0)),
#'                       stage1 = 1:2, w1 = 0.7, cac_plan = 0.8)
#' alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
#' set.seed(1)
#' effect <- rnorm(20, 0, 0.3)
#' data_all <- expand.grid(cluster = 1:20, period = 1:4)
#' data_all$m <- 40
#' data_all$y <- rbinom(80, 40, plogis(-1 + effect[data_all$cluster]))
#' interim <- interim_statistic(plan, data_all[data_all$period <= 2, ], alloc)
#' final <- final_statistic(interim, data_all)
#' effect_estimates(final)
#' @export
effect_estimates <- function(object, check = TRUE) {
  final <- inherits(object, "final_statistic")
  if (!final && !inherits(object, "interim_statistic")) {
    stop("object must be an interim_statistic or final_statistic object")
  }
  interim <- if (final) object$interim else object
  plan <- interim$plan
  nuis <- interim$nuis
  fam <- plan$fam
  fin <- if (final) object[c("design", "data")] else NULL

  stat <- function(k) function(b) .stage_stats(interim, fin, b, k)[[k]]
  z_crit <- stats::qnorm(1 - plan$alpha / 2)
  rows <- if (final) {
    list(adaptive = list("Z", plan$c2_ci), stage1 = list("z1", z_crit),
         stage2 = list("z2", z_crit), unadjusted = list("z_all", z_crit))
  } else {
    list(adaptive = list("z1", plan$c1), unadjusted = list("z1", z_crit))
  }
  scale <- if (plan$family == "gaussian") sqrt(nuis$A + nuis$sigma2) else 1

  fits <- do.call(rbind, lapply(rows, function(r) {
    .invert_statistic(stat(r[[1]]), r[[2]], scale)
  }))
  res <- data.frame(
    analysis = names(rows),
    critical_value = vapply(rows, function(r) r[[2]], numeric(1)),
    fits,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (anyNA(res$estimate)) {
    warning("A point estimate could not be found: ",
            paste(res$analysis[is.na(res$estimate)], collapse = ", "),
            call. = FALSE)
  }
  unbounded <- is.infinite(res$lower) | is.infinite(res$upper)
  if (any(unbounded & is.finite(res$critical_value))) {
    warning("An interval limit could not be bounded: ",
            paste(res$analysis[unbounded & is.finite(res$critical_value)],
                  collapse = ", "),
            call. = FALSE)
  }

  if (fam$has_ratio) {
    res$ratio <- exp(res$estimate)
    res$ratio_lower <- exp(res$lower)
    res$ratio_upper <- exp(res$upper)
  } else {
    res$ratio <- res$ratio_lower <- res$ratio_upper <- NA_real_
  }

  source_obj <- if (final) object else interim
  omega <- rowsum(ifelse(source_obj$data$obs, source_obj$data$M, 0),
                  source_obj$design$strata)
  omega <- omega[rownames(omega) %in% nuis$strata$stratum, , drop = FALSE]
  omega <- omega / sum(omega)

  difference <- function(b) {
    vapply(b, function(x) {
      if (!is.finite(x)) return(NA_real_)
      pars <- .linear_pars(plan, nuis, x, source_obj$data,
                           source_obj$design$Tn)
      base <- outer(pars$eta[rownames(omega)], pars$gamma, "+")
      used <- omega > 0
      sum(omega[used] * (fam$marginal(base[used] + x, nuis$A, nuis$gh) -
                           fam$marginal(base[used], nuis$A, nuis$gh)))
    }, numeric(1))
  }
  h <- 1e-3 * scale
  res$difference <- difference(res$estimate)
  res$difference_lower <- difference(res$lower)
  res$difference_upper <- difference(res$upper)
  res$difference_se <- res$se *
    (difference(res$estimate + h) - difference(res$estimate - h)) / (2 * h)

  if (check) {
    a <- res[res$analysis == "adaptive", ]
    if (is.finite(a$lower) && is.finite(a$upper)) {
      width <- a$upper - a$lower
      g <- seq(a$lower - width / 2, a$upper + width / 2, length.out = 21)
      values <- vapply(g, stat(rows$adaptive[[1]]), numeric(1))
      if (any(diff(values) >= 0)) {
        warning("Primary statistic is not monotone in the treatment effect: ",
                "inspect before reporting.", call. = FALSE)
      }
    }
  }

  if (!isTRUE(all.equal(plan$c2, plan$c2_ci)) && final) {
    message(sprintf(
      "Repeated CI uses the non-binding critical value %.3f (test boundary: %.3f).",
      plan$c2_ci, plan$c2
    ))
  }

  attr(res, "family") <- plan$family
  attr(res, "labels") <- fam$labels
  attr(res, "stage") <- if (final) "final" else "interim"
  class(res) <- c("effect_estimates", "data.frame")
  res
}

#' @export
print.effect_estimates <- function(x, digits = 3, ...) {
  labels <- attr(x, "labels")
  stage <- attr(x, "stage")
  cat(sprintf("Treatment effect estimates (%s analysis, %s)\n",
              stage, attr(x, "family")))

  format_row <- function(est, lo, up) {
    sprintf("%s (%s, %s)", format(round(est, digits), nsmall = digits),
            format(round(lo, digits), nsmall = digits),
            format(round(up, digits), nsmall = digits))
  }
  df <- data.frame(analysis = x$analysis, stringsAsFactors = FALSE)
  df[[labels[["link"]]]] <- format_row(x$estimate, x$lower, x$upper)
  df$se <- format(round(x$se, digits), nsmall = digits)
  if (!is.na(labels[["ratio"]])) {
    df[[labels[["ratio"]]]] <- format_row(x$ratio, x$ratio_lower, x$ratio_upper)
    df[[labels[["difference"]]]] <- format_row(x$difference, x$difference_lower,
                                               x$difference_upper)
  }
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}
