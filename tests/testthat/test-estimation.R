schedule_xo <- rbind(c(0, 1, 0, 1), c(1, 0, 1, 0))

simulate_trial <- function(family, k = 12, n_strata = 3, m = 40, beta = 0.3,
                           b0 = -1, A = 0.2, cac = 0.7, seed = 1) {
  set.seed(seed)
  N <- n_strata * 2 * k
  alloc <- data.frame(
    cluster = seq_len(N),
    country = paste0("C", rep(seq_len(n_strata), each = 2 * k)),
    sequence = rep(rep(1:2, each = k), n_strata)
  )
  u <- matrix(0, N, 4)
  u[, 1] <- rnorm(N, 0, sqrt(A))
  for (t in 2:4) u[, t] <- cac * u[, t - 1] + rnorm(N, 0, sqrt(A * (1 - cac^2)))
  eta <- b0 + beta * schedule_xo[alloc$sequence, ] + u
  cp <- data.frame(
    cluster = rep(seq_len(N), 4), period = rep(1:4, each = N),
    m = rpois(4 * N, m) + 5L, eta = c(eta)
  )
  cp$y <- switch(
    family,
    binomial = rbinom(nrow(cp), cp$m, plogis(cp$eta)),
    poisson = rpois(nrow(cp), cp$m * exp(cp$eta)),
    gaussian = rnorm(nrow(cp), cp$eta, 1 / sqrt(cp$m)) * cp$m
  )
  if (family == "gaussian") cp$s2 <- rchisq(nrow(cp), cp$m - 1) / (cp$m - 1)
  cp$eta <- NULL
  list(alloc = alloc, cp = cp)
}

make_plan <- function(family, ...) {
  analysis_plan(family = family, schedule = schedule_xo, stage1 = 1:2,
                strata = "country", w1 = 0.7, cac_plan = 0.8, ...)
}

run_trial <- function(family, ...) {
  tr <- simulate_trial(family, ...)
  plan <- make_plan(family)
  int <- interim_statistic(plan, tr$cp[tr$cp$period <= 2, ], tr$alloc)
  list(tr = tr, plan = plan, int = int, fin = final_statistic(int, tr$cp))
}

test_that("stage_statistics at beta = 0 reproduces the test statistics", {
  for (family in c("binomial", "poisson", "gaussian")) {
    r <- run_trial(family)
    s <- stage_statistics(r$fin, 0)
    expect_equal(s$z1, r$int$z1, tolerance = 1e-10)
    expect_equal(s$z2, r$fin$z2, tolerance = 1e-10)
    expect_equal(s$Z, r$fin$Z, tolerance = 1e-10)
    expect_equal(r$fin$Z, r$plan$w1 * r$fin$z1 + r$plan$w2 * r$fin$z2)
  }
})

test_that("statistics decrease in beta and are zero at the estimate", {
  for (family in c("binomial", "poisson", "gaussian")) {
    r <- run_trial(family)
    z <- vapply(c(-0.5, 0, 0.5, 1), function(b) stage_statistics(r$fin, b)$Z, 0)
    expect_true(all(diff(z) < 0))
    est <- effect_estimates(r$fin)
    adaptive <- est[est$analysis == "adaptive", ]
    expect_equal(stage_statistics(r$fin, adaptive$estimate)$Z, 0, tolerance = 1e-5)
    expect_equal(stage_statistics(r$fin, adaptive$lower)$Z, r$plan$c2,
                 tolerance = 1e-5)
    expect_equal(stage_statistics(r$fin, adaptive$upper)$Z, -r$plan$c2,
                 tolerance = 1e-5)
  }
})

test_that("the adaptive estimate is close to the truth in a large trial", {
  for (family in c("binomial", "poisson", "gaussian")) {
    r <- run_trial(family, k = 80, seed = 2)
    est <- effect_estimates(r$fin, check = FALSE)
    adaptive <- est[est$analysis == "adaptive", ]
    expect_lt(abs(adaptive$estimate - 0.3), 3 * adaptive$se)
  }
})

test_that("effect scales follow the family", {
  est_b <- effect_estimates(run_trial("binomial")$fin)
  expect_equal(est_b$ratio, exp(est_b$estimate))
  est_g <- effect_estimates(run_trial("gaussian")$fin)
  expect_true(all(is.na(est_g$ratio)))
  expect_equal(est_g$difference, est_g$estimate, tolerance = 1e-10)
})

test_that("the marginal difference weights each stratum separately", {
  r <- run_trial("binomial", n_strata = 4)
  nuis <- r$int$nuis
  eta <- nuis$eta
  by_stratum <- r$plan$fam$marginal(eta, nuis$A, nuis$gh)
  expect_length(by_stratum, 4)
  expect_equal(unname(by_stratum[1]),
               r$plan$fam$marginal(eta[[1]], nuis$A, nuis$gh))
})

test_that("a trial stopped at the interim gets interim rows only", {
  r <- run_trial("binomial")
  est <- effect_estimates(r$int)
  expect_equal(est$analysis, c("adaptive", "unadjusted"))
  expect_equal(est$critical_value[1], r$plan$c1)
})

test_that("the gaussian posterior matches a fine grid", {
  set.seed(3)
  n <- 6
  M <- matrix(rpois(4 * n, 30), n, 4)
  obs <- M > 0
  ETA <- matrix(rnorm(4 * n, 0, 0.2), n, 4)
  Y <- M * (ETA + rnorm(4 * n, 0, 0.3))
  fam <- list(linkinv = identity,
              loglik = function(H, y, m) -0.5 * m * (y / m - H)^2 / 0.9)
  exact <- adaptiveCRT:::.posterior_mean_gaussian(Y, M, ETA, obs, 0.1, 0.7, 0.9)
  grid <- adaptiveCRT:::.posterior_grid(0.1, 0.7, list(span = 7, grid_min = 1500L))
  expect_equal(exact, adaptiveCRT:::.posterior_mean_grid(fam, Y, M, ETA, obs, grid),
               tolerance = 1e-6)
})

test_that("individual-level data give the same result as cluster-period data", {
  set.seed(4)
  alloc <- data.frame(cluster = 1:12, sequence = rep(1:2, 6))
  ind <- expand.grid(cluster = 1:12, period = 1:4, id = 1:20)
  ind$outcome <- rbinom(nrow(ind), 1, 0.3)
  ind$outcome[sample(nrow(ind), 30)] <- NA
  obs <- ind[!is.na(ind$outcome), ]
  cp <- aggregate(outcome ~ cluster + period, obs, function(v) c(sum(v), length(v)))
  cp <- data.frame(cluster = cp$cluster, period = cp$period,
                   y = cp$outcome[, 1], m = cp$outcome[, 2])
  plan <- analysis_plan(family = "binomial", schedule = schedule_xo, stage1 = 1:2,
                        w1 = 0.7, cac_plan = 0.8)
  a <- interim_statistic(plan, ind[ind$period <= 2, ], alloc)
  b <- interim_statistic(plan, cp[cp$period <= 2, ], alloc)
  expect_equal(a$z1, b$z1)
  expect_equal(a$nuis$A, b$nuis$A)
})

test_that("a single stage 1 period uses the planned autocorrelation", {
  set.seed(5)
  sched <- rbind(c(0, 0), c(1, 1))
  alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
  cp <- expand.grid(cluster = 1:20, period = 1:2)
  cp$m <- 30
  cp$y <- rbinom(40, 30, 0.3)
  plan <- analysis_plan(family = "binomial", schedule = sched, stage1 = 1,
                        w1 = 0.7, cac_plan = 0.65)
  int <- interim_statistic(plan, cp[cp$period == 1, ], alloc)
  expect_equal(int$nuis$cac, 0.65)
  expect_true(is.finite(final_statistic(int, cp)$Z))
})

test_that("clusters added at stage 2 have no stage 1 contribution", {
  set.seed(6)
  sched <- rbind(c(0, 0), c(1, 1), c(NA, 0), c(NA, 1))
  alloc <- data.frame(cluster = 1:20, sequence = rep(1:2, each = 10))
  alloc2 <- rbind(alloc, data.frame(cluster = 21:24, sequence = rep(3:4, each = 2)))
  cp <- rbind(
    data.frame(cluster = 1:20, period = 1, m = 30, y = rbinom(20, 30, 0.3)),
    data.frame(cluster = 1:24, period = 2, m = 30, y = rbinom(24, 30, 0.3))
  )
  plan <- analysis_plan(family = "binomial", schedule = sched, stage1 = 1,
                        w1 = 0.7, cac_plan = 0.8)
  int <- interim_statistic(plan, cp[cp$period == 1, ], alloc)
  expect_error(final_statistic(int, cp), "not in the allocation")
  fin <- final_statistic(int, cp, alloc = alloc2)
  expect_length(fin$U2, 24)
  expect_true(is.finite(fin$Z))

  bad <- alloc2
  bad$sequence[21:22] <- 1L
  expect_error(final_statistic(int, cp, alloc = bad), "NA in the schedule")
})

test_that("the grid is refined for very informative cluster-periods", {
  fam <- adaptiveCRT:::.estimation_family("poisson")
  plan <- list(fam = fam, control = adaptiveCRT:::.estimation_control(list()))
  nuis <- list(A = 0.09, cac = 0.7,
               grid = adaptiveCRT:::.posterior_grid(0.09, 0.7, plan$control))
  small <- list(Y = matrix(10, 2, 2), M = matrix(50, 2, 2),
                obs = matrix(TRUE, 2, 2))
  large <- list(Y = matrix(5000, 2, 2), M = matrix(20000, 2, 2),
                obs = matrix(TRUE, 2, 2))
  expect_length(adaptiveCRT:::.resolve_grid(plan, nuis, small)$b, 80)
  expect_gt(length(adaptiveCRT:::.resolve_grid(plan, nuis, large)$b), 80)
})

test_that("bad input is rejected", {
  tr <- simulate_trial("binomial", k = 4, n_strata = 1)
  plan <- make_plan("binomial")
  expect_error(interim_statistic(plan, tr$cp, tr$alloc[-1, ]), "not in the allocation")
  bad <- tr$cp
  bad$y[1] <- bad$m[1] + 1
  expect_error(interim_statistic(plan, bad[bad$period <= 2, ], tr$alloc), "between 0 and m")
  expect_error(analysis_plan(family = "weibull", schedule = schedule_xo,
                             stage1 = 1:2, w1 = 0.7, cac_plan = 0.8), "family must be")
  expect_error(analysis_plan(family = "binomial", schedule = schedule_xo,
                             stage1 = 1:4, w1 = 0.7, cac_plan = 0.8), "after stage 1")

  gaussian_plan <- make_plan("gaussian")
  gd <- simulate_trial("gaussian", k = 4, n_strata = 1)
  gd$cp$s2 <- NULL
  expect_error(interim_statistic(gaussian_plan, gd$cp[gd$cp$period <= 2, ], gd$alloc),
               "sigma2")
})

simulate_wedge <- function(family = "binomial", beta = 0.3, trend = seq(0, 1, length.out = 6),
                           k = 3, seed = 8) {
  set.seed(seed)
  sched <- t(vapply(1:5, function(s) as.numeric(1:6 > s), numeric(6)))
  alloc <- data.frame(cluster = seq_len(5 * k), sequence = rep(1:5, each = k))
  u <- rnorm(5 * k, 0, 0.3)
  cp <- expand.grid(cluster = alloc$cluster, period = 1:6)
  cp$m <- 40
  eta <- -1 + trend[cp$period] + beta * sched[cbind(alloc$sequence[cp$cluster], cp$period)] +
    u[cp$cluster]
  cp$y <- switch(
    family,
    binomial = rbinom(nrow(cp), cp$m, plogis(eta)),
    poisson = rpois(nrow(cp), cp$m * exp(eta)),
    gaussian = rnorm(nrow(cp), eta, 1 / sqrt(cp$m)) * cp$m
  )
  if (family == "gaussian") cp$s2 <- 1
  list(schedule = sched, alloc = alloc, cp = cp)
}

test_that("period effects default to the design's balance", {
  expect_equal(make_plan("binomial")$period_effects, "none")
  wedge <- simulate_wedge()
  plan <- analysis_plan(family = "binomial", schedule = wedge$schedule, stage1 = 1:3,
                        w1 = 0.7, cac_plan = 0.8)
  expect_equal(plan$period_effects, "fixed")
  override <- analysis_plan(family = "binomial", schedule = wedge$schedule, stage1 = 1:3,
                            period_effects = "none", w1 = 0.7, cac_plan = 0.8)
  expect_equal(override$period_effects, "none")
})

test_that("period effects are recovered and match the statistics at beta = 0", {
  for (family in c("binomial", "poisson", "gaussian")) {
    wedge <- simulate_wedge(family)
    plan <- analysis_plan(family = family, schedule = wedge$schedule, stage1 = 1:3,
                          w1 = 0.7, cac_plan = 0.8)
    int <- interim_statistic(plan, wedge$cp[wedge$cp$period <= 3, ], wedge$alloc)
    fin <- final_statistic(int, wedge$cp)
    s <- stage_statistics(fin, 0)
    expect_equal(s$z1, int$z1, tolerance = 1e-10)
    expect_equal(s$Z, fin$Z, tolerance = 1e-10)

    pars <- adaptiveCRT:::.linear_pars(plan, int$nuis, 0.3, fin$data, fin$design$Tn)
    expect_equal(pars$gamma[1], 0)
    expect_lt(max(abs(pars$gamma - seq(0, 1, length.out = 6))), 0.35)
  }
})

test_that("period effects correct a stepped-wedge estimate under a strong trend", {
  wedge <- simulate_wedge("gaussian", k = 10, trend = c(0, 1, 2, 3, 4, 5))
  fit <- function(pe) {
    plan <- analysis_plan(family = "gaussian", schedule = wedge$schedule, stage1 = 1:3,
                          period_effects = pe, w1 = 0.7, cac_plan = 0.8)
    int <- suppressWarnings(interim_statistic(plan, wedge$cp[wedge$cp$period <= 3, ],
                                              wedge$alloc))
    effect_estimates(final_statistic(int, wedge$cp), check = FALSE)
  }
  fixed <- fit("fixed")
  expect_lt(abs(fixed$estimate[fixed$analysis == "adaptive"] - 0.3), 0.15)
})

test_that("one cluster per sequence can use stage 1 pattern for variance estimation", {
  wedge <- simulate_wedge(k = 1)
  wedge$schedule <- t(vapply(1:5, function(s) as.numeric(1:6 > s), numeric(6)))
  plan <- analysis_plan(family = "binomial", schedule = wedge$schedule, stage1 = 1:3,
                        w1 = 0.7, cac_plan = 0.8)
  int <- interim_statistic(plan, wedge$cp[wedge$cp$period <= 3, ], wedge$alloc)
  expect_equal(int$nuis$raw$df, 2)
})

test_that("a changed stage 2 schedule is accepted only if stage 1 is unchanged", {
  wedge <- simulate_wedge()
  plan <- analysis_plan(family = "binomial", schedule = wedge$schedule, stage1 = 1:3,
                        w1 = 0.7, cac_plan = 0.8)
  int <- interim_statistic(plan, wedge$cp[wedge$cp$period <= 3, ], wedge$alloc)
  shorter <- wedge$schedule[, 1:5]
  fin <- final_statistic(int, wedge$cp[wedge$cp$period <= 5, ], schedule = shorter)
  expect_equal(fin$design$Tn, 5)
  expect_true(is.finite(fin$Z))

  changed <- wedge$schedule
  changed[1, 2] <- 0
  expect_error(final_statistic(int, wedge$cp, schedule = changed), "stage 1 treatment")
})

test_that("stratification is optional and may use several factors", {
  tr <- simulate_trial("binomial", n_strata = 2, k = 8)
  cp1 <- tr$cp[tr$cp$period <= 2, ]

  plan_none <- analysis_plan(family = "binomial", schedule = schedule_xo, stage1 = 1:2,
                             w1 = 0.7, cac_plan = 0.8)
  single <- transform(tr$alloc, country = "only")
  plan_single <- make_plan("binomial")
  expect_equal(interim_statistic(plan_none, cp1, tr$alloc)$z1,
               interim_statistic(plan_single, cp1, single)$z1)

  alloc <- tr$alloc
  alloc$size <- rep(rep(c("small", "large"), each = 4), 4)
  plan_two <- analysis_plan(family = "binomial", schedule = schedule_xo, stage1 = 1:2,
                            strata = c("country", "size"), w1 = 0.7, cac_plan = 0.8)
  int <- interim_statistic(plan_two, cp1, alloc)
  expect_setequal(int$nuis$strata$stratum,
                  c("C1:small", "C1:large", "C2:small", "C2:large"))
  expect_true(is.finite(int$z1))
  expect_s3_class(effect_estimates(int), "effect_estimates")
})

test_that("strata too small for variance estimation fall back with a warning", {
  tr <- simulate_trial("binomial", n_strata = 1, k = 8)
  alloc <- tr$alloc
  alloc$pair <- rep(1:8, 2)
  plan <- analysis_plan(family = "binomial", schedule = schedule_xo, stage1 = 1:2,
                        strata = "pair", w1 = 0.7, cac_plan = 0.8)
  expect_warning(int <- interim_statistic(plan, tr$cp[tr$cp$period <= 2, ], alloc),
                 "too small")
  expect_equal(int$nuis$raw$df, 14)
})

test_that("a stratum with no contrast is flagged", {
  tr <- simulate_trial("binomial", n_strata = 1, k = 8)
  alloc <- tr$alloc
  alloc$site <- c(rep("A", 15), "B")
  plan <- analysis_plan(family = "binomial", schedule = schedule_xo, stage1 = 1:2,
                        strata = "site", w1 = 0.7, cac_plan = 0.8)
  expect_warning(interim_statistic(plan, tr$cp[tr$cp$period <= 2, ], alloc),
                 "No treatment contrast")
})

test_that("interim_decision needs a linked design", {
  r <- run_trial("binomial", k = 4, n_strata = 2)
  expect_error(interim_decision(r$int), "not linked")
})

test_that("a design object supplies the plan", {
  skip_on_cran()
  design <- parallel_crt(icc = 0.05, cac = 0.8, delta = -0.1, k1 = 12, m1 = 30,
                         k2 = 0:2, m2 = c(20, 30), rho = 30,
                         family = "binomial", baseline = 0.3)
  planned <- suppressWarnings(
    adaptive_analysis(design, target_power = 0.6, verbose = FALSE, tol = 0.02,
                      method = "lambda")
  )
  plan <- analysis_plan(planned)
  params <- planned$raw$results$params
  expect_equal(plan$family, "binomial")
  expect_equal(plan$w1, params$w1_ref)
  expect_equal(plan$c1, params$efficacy_boundary)
  expect_equal(plan$c2, params$c2)
  expect_equal(plan$cac_plan, 0.8)
  expect_equal(nrow(plan$schedule), 4)
  expect_error(analysis_plan(planned, family = "poisson"), "differs")

  set.seed(7)
  alloc <- data.frame(cluster = 1:24, sequence = rep(1:2, each = 12))
  cp <- data.frame(cluster = 1:24, period = 1, m = 30, y = rbinom(24, 30, 0.3))
  int <- interim_statistic(plan, cp, alloc)
  dec <- interim_decision(int, verbose = FALSE)
  expect_s3_class(dec, "interim_result")
  expect_equal(dec$z1, int$z1)
  expect_warning(interim_statistic(plan, cp[-1, ], alloc[-1, ]), "planned 24")
})
