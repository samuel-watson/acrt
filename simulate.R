#' Simulate adaptive CRT with individual-level data and lmer estimation
#'
#' @param results Output from adaptive_analysis (single design, not explored)
#' @param n_sims Number of simulations
#' @param hypothesis "H0" for type I error, "H1" for power
#' @param seed Random seed
#' @param verbose Print progress
#' @return List with empirical and theoretical results
simulate_adaptive_trial <- function(results,
                                    n_sims = 1000,
                                    hypothesis = c("H1", "H0"),
                                    seed = 12345,
                                    verbose = TRUE) {

  hypothesis <- match.arg(hypothesis)
  set.seed(seed)

  # === Extract design parameters ===
  raw <- extract_raw_results(results, 1)

  if (is.null(raw)) {
    stop("Could not extract results. Use a single-design adaptive_analysis result.")
  }

  design <- results$design
  spec <- design$spec
  fixed_params <- design$fixed_params
  opt_designs <- raw$quadrature$optimal_designs

  # True treatment effect
  delta_true <- if (hypothesis == "H1") fixed_params$delta else 0

  # Variance components
  icc <- fixed_params$icc
  sigma2_total <- 1  # Assume standardised
  sigma2_b <- icc * sigma2_total  # Between-cluster variance
  sigma2_w <- (1 - icc) * sigma2_total  # Within-cluster variance

  if (verbose) {
    cat("=== Simulation Setup ===\n")
    cat(sprintf("Hypothesis: %s\n", hypothesis))
    cat(sprintf("True delta: %.3f\n", delta_true))
    cat(sprintf("ICC: %.3f (sigma2_b=%.3f, sigma2_w=%.3f)\n", icc, sigma2_b, sigma2_w))
    cat(sprintf("Design type: %s\n", spec$design_type))
    cat(sprintf("Simulations: %d\n\n", n_sims))
  }

  # === Build decision rule lookup ===
  decision_lookup <- opt_designs
  design_vars <- setdiff(names(opt_designs),
                         c("z1", "efficacy_stop", "futility_stop", "continue",
                           "cp", "criterion", "cost"))

  # === Get model structures via model builder ===
  if (verbose) cat("Building model structures...\n")

  design_grid <- raw$models$design_grid
  stage1_params <- design$stage1_grid[1, , drop = FALSE]

  # Build base fixed params
  base_fixed_params <- c(as.list(fixed_params), as.list(stage1_params))
  base_fixed_params$cache <- new.env()

  # Get stage 1 design info from first model output
  model_output_1 <- raw$models$list[[1]]
  k1 <- model_output_1$k1
  m1 <- model_output_1$m1
  n_arms <- 2

  if (verbose) {
    cat(sprintf("Stage 1: k1=%d clusters/arm, m1=%d individuals/cluster\n", k1, m1))
  }

  # Pre-compute stage 2 model outputs for weights
  # Pre-compute stage 2 weights lookup - keyed by (m2, k2)
  # These are fixed based on design, not data
  stage2_weights <- list()

  for (i in seq_len(nrow(design_grid))) {
    s2_params <- design_grid[i, , drop = FALSE]
    base_fixed_params$cache$mod <- NULL
    base_fixed_params$cache$cache_key <- NULL

    tryCatch({
      mod_out <- spec$model_builder(s2_params, base_fixed_params)

      # Create key from design parameters
      m2_val <- s2_params$m2 %||% m1
      k2_val <- s2_params$k2 %||% 0
      key <- sprintf("m2=%g_k2=%g", m2_val, k2_val)

      stage2_weights[[key]] <- list(
        w1 = mod_out$w1,
        w2 = mod_out$w2,
        m2 = m2_val,
        k2 = k2_val
      )
    }, error = function(e) {
      if (verbose) cat(sprintf("  Model %d error: %s\n", i, e$message))
    })
  }

  if (verbose) {
    cat(sprintf("Pre-computed weights for %d stage 2 designs:\n", length(stage2_weights)))
    for (key in names(stage2_weights)) {
      w <- stage2_weights[[key]]
      cat(sprintf("  %s: w1=%.4f, w2=%.4f (w1^2+w2^2=%.4f)\n",
                  key, w$w1, w$w2, w$w1^2 + w$w2^2))
    }
    cat("\n")
  }

  # === Simulation function ===
  simulate_one_trial <- function(sim_id) {

    # =========================================
    # STAGE 1: Generate individual-level data
    # =========================================

    # Create cluster IDs
    # Control arm: clusters 1 to k1
    # Treatment arm: clusters (k1+1) to (2*k1)
    cluster_ids_ctrl <- 1:k1
    cluster_ids_trt <- (k1 + 1):(2 * k1)

    # Generate cluster random effects
    u_ctrl <- rnorm(k1, 0, sqrt(sigma2_b))
    u_trt <- rnorm(k1, 0, sqrt(sigma2_b))

    # Store for stage 2 continuity
    cluster_effects <- c(u_ctrl, u_trt)
    names(cluster_effects) <- c(cluster_ids_ctrl, cluster_ids_trt)

    # Generate individual-level data for control arm
    df_ctrl <- data.frame(
      id = 1:(k1 * m1),
      cl = rep(cluster_ids_ctrl, each = m1),
      trt = 0
    )
    df_ctrl$y <- 0 + # intercept
      cluster_effects[as.character(df_ctrl$cl)] +
      rnorm(nrow(df_ctrl), 0, sqrt(sigma2_w))

    # Generate individual-level data for treatment arm
    df_trt <- data.frame(
      id = (k1 * m1 + 1):(2 * k1 * m1),
      cl = rep(cluster_ids_trt, each = m1),
      trt = 1
    )
    df_trt$y <- delta_true + # intercept + treatment effect
      cluster_effects[as.character(df_trt$cl)] +
      rnorm(nrow(df_trt), 0, sqrt(sigma2_w))

    # Combine stage 1 data
    df_s1 <- rbind(df_ctrl, df_trt)
    df_s1$cl <- factor(df_s1$cl)
    df_s1$stage <- 1

    # =========================================
    # STAGE 1: Fit mixed model with lmer
    # =========================================

    # =========================================
    # STAGE 1: Fit mixed model with lmerTest
    # =========================================

    fit_s1 <- tryCatch({
      lmerTest::lmer(y ~ trt + (1|cl), data = df_s1, REML = TRUE)
    }, error = function(e) {
      NULL
    })

    if (is.null(fit_s1)) {
      return(list(reject = NA, stage = NA, z1 = NA, error = "Stage 1 lmer failed"))
    }

    # Extract treatment effect estimate with Satterthwaite df
    coef_s1 <- summary(fit_s1)$coefficients
    theta1_hat <- coef_s1["trt", "Estimate"]
    se1_hat <- coef_s1["trt", "Std. Error"]
    df1 <- coef_s1["trt", "df"]
    t1 <- theta1_hat / se1_hat

    # For combination test, we still need a z-scale statistic
    # Convert t to z-scale using probability integral transform
    p1 <- pt(t1, df = df1)  # CDF of t
    z1 <- qnorm(p1)         # Convert to z-scale

    # =========================================
    # INTERIM DECISION
    # =========================================

    # Find closest z1 in decision lookup
    idx <- which.min(abs(decision_lookup$z1 - z1))
    decision <- decision_lookup[idx, ]

    # Check for efficacy stop
    # Check for efficacy stop - use WEIGHTED statistic
    if (decision$efficacy_stop) {
      # Get w1 for this decision (use first available as reference)
      w1_eff <- stage2_weights[[names(stage2_weights)[1]]]$w1

      # Weighted stage 1 statistic
      z1_weighted <- w1_eff * z1

      # Reject only if weighted statistic exceeds critical value
      reject_s1 <- abs(z1_weighted) > 1.96

      return(list(
        reject = reject_s1,
        stage = 1,
        z1 = z1,
        z1_weighted = z1_weighted,
        z2 = NA,
        z_combo = NA,
        w1 = w1_eff,
        w2 = NA,
        decision = "efficacy",
        theta1_hat = theta1_hat,
        se1_hat = se1_hat,
        t1 = t1,
        df1 = df1
      ))
    }

    # Check for futility stop
    if (decision$futility_stop) {
      return(list(
        reject = FALSE,
        stage = 1,
        z1 = z1,
        z2 = NA,
        z_combo = NA,
        w1 = NA,
        w2 = NA,
        decision = "futility",
        theta1_hat = theta1_hat,
        se1_hat = se1_hat,
        t1 = t1,
        df1 = df1
      ))
    }

    # =========================================
    # STAGE 2: Continue with adapted design
    # =========================================

    # Get adapted design parameters from decision rule
    m2 <- decision$m2
    k2 <- decision$k2 %||% 0

    # Handle NA values
    m2_actual <- if (is.na(m2)) m1 else m2
    k2_actual <- if (is.na(k2)) 0 else k2

    # Look up PRE-COMPUTED weights based on design choice (not data!)
    weight_key <- sprintf("m2=%g_k2=%g", m2_actual, k2_actual)

    if (!weight_key %in% names(stage2_weights)) {
      # Find closest match
      available_keys <- names(stage2_weights)
      if (length(available_keys) == 0) {
        return(list(reject = NA, stage = 2, z1 = z1, z2 = NA,
                    z_combo = NA, decision = "continue",
                    error = "No stage 2 weights available"))
      }
      weight_key <- available_keys[1]
      if (sim_id == 1) {
        cat(sprintf("Warning: No exact weight match for m2=%g, k2=%g. Using %s\n",
                    m2_actual, k2_actual, weight_key))
      }
    }

    # Get fixed weights
    s2_weights <- stage2_weights[[weight_key]]
    w1 <- s2_weights$w1
    w2 <- s2_weights$w2

    # Total clusters in stage 2: original k1 per arm + k2 new per arm
    total_k_s2_per_arm <- k1 + k2_actual

    # =========================================
    # STAGE 2: Generate data
    # =========================================

    # New cluster IDs for stage 2 (if any)
    if (k2_actual > 0) {
      new_cluster_ids_ctrl <- (2 * k1 + 1):(2 * k1 + k2_actual)
      new_cluster_ids_trt <- (2 * k1 + k2_actual + 1):(2 * k1 + 2 * k2_actual)

      # Generate new cluster effects
      u_new_ctrl <- rnorm(k2_actual, 0, sqrt(sigma2_b))
      u_new_trt <- rnorm(k2_actual, 0, sqrt(sigma2_b))

      names(u_new_ctrl) <- new_cluster_ids_ctrl
      names(u_new_trt) <- new_cluster_ids_trt

      cluster_effects <- c(cluster_effects, u_new_ctrl, u_new_trt)
    } else {
      new_cluster_ids_ctrl <- integer(0)
      new_cluster_ids_trt <- integer(0)
    }

    # All stage 2 clusters
    all_cluster_ids_ctrl_s2 <- c(cluster_ids_ctrl, new_cluster_ids_ctrl)
    all_cluster_ids_trt_s2 <- c(cluster_ids_trt, new_cluster_ids_trt)

    n_clusters_s2_per_arm <- length(all_cluster_ids_ctrl_s2)

    # Generate stage 2 control data
    df_s2_ctrl <- data.frame(
      id = max(df_s1$id) + 1:(n_clusters_s2_per_arm * m2_actual),
      cl = rep(all_cluster_ids_ctrl_s2, each = m2_actual),
      trt = 0
    )
    df_s2_ctrl$y <- 0 +
      cluster_effects[as.character(df_s2_ctrl$cl)] +
      rnorm(nrow(df_s2_ctrl), 0, sqrt(sigma2_w))

    # Generate stage 2 treatment data
    df_s2_trt <- data.frame(
      id = max(df_s2_ctrl$id) + 1:(n_clusters_s2_per_arm * m2_actual),
      cl = rep(all_cluster_ids_trt_s2, each = m2_actual),
      trt = 1
    )
    df_s2_trt$y <- delta_true +
      cluster_effects[as.character(df_s2_trt$cl)] +
      rnorm(nrow(df_s2_trt), 0, sqrt(sigma2_w))

    # Combine stage 2 data
    df_s2 <- rbind(df_s2_ctrl, df_s2_trt)
    df_s2$cl <- factor(df_s2$cl)
    df_s2$stage <- 2

    # =========================================
    # FINAL TEST using full data
    # =========================================

    # Combine all data
    df_full <- rbind(df_s1, df_s2)
    df_full$cl <- factor(df_full$cl)

    fit_full <- tryCatch({
      lmerTest::lmer(y ~ trt + (1|cl), data = df_full, REML = TRUE)
    }, error = function(e) {
      NULL
    })

    if (is.null(fit_full)) {
      return(list(reject = NA, stage = 2, z1 = z1, z2 = NA,
                  z_combo = NA, w1 = w1, w2 = w2,
                  decision = "continue", error = "Full model failed"))
    }

    # Final test based on full data
    coef_full <- summary(fit_full)$coefficients
    theta_full <- coef_full["trt", "Estimate"]
    se_full <- coef_full["trt", "Std. Error"]
    df_full_sat <- coef_full["trt", "df"]
    t_full <- theta_full / se_full

    # Two-sided p-value with t-distribution
    p_full <- 2 * pt(abs(t_full), df = df_full_sat, lower.tail = FALSE)
    reject <- p_full < 0.05

    # For diagnostics: convert to z-scale
    z_combo <- sign(t_full) * qnorm(1 - p_full/2)

    # Back-calculate implied z2|1 for diagnostics only
    z2 <- (z_combo - w1 * z1) / w2

    list(
      reject = reject,
      stage = 2,
      z1 = z1,
      z2 = z2,
      z_combo = z_combo,
      w1 = w1,
      w2 = w2,
      decision = "continue",
      theta1_hat = theta1_hat,
      theta_full = theta_full,
      se1_hat = se1_hat,
      se_full = se_full,
      t1 = t1,
      t_full = t_full,
      df1 = df1,
      df_full = df_full_sat,
      p_full = p_full,
      m2 = m2_actual,
      k2 = k2_actual
    )
  }

  # === Run simulations ===
  if (verbose) cat("Running simulations...\n")

  sim_results <- vector("list", n_sims)
  errors <- 0

  for (i in seq_len(n_sims)) {
    if (verbose && i %% 200 == 0) {
      cat(sprintf("  %d / %d (errors: %d)\n", i, n_sims, errors))
    }

    result <- tryCatch({
      simulate_one_trial(i)
    }, error = function(e) {
      list(reject = NA, stage = NA, z1 = NA, error = e$message)
    })

    if (!is.null(result$error)) errors <- errors + 1
    sim_results[[i]] <- result
  }

  # === Summarise results ===
  valid_results <- sim_results[!sapply(sim_results, function(x) is.na(x$reject))]
  n_valid <- length(valid_results)

  if (n_valid == 0) {
    # Print some error messages for debugging
    error_msgs <- sapply(sim_results[1:min(5, n_sims)], function(x) x$error %||% "unknown")
    cat("First few errors:\n")
    print(error_msgs)
    stop("All simulations failed")
  }

  rejections <- sapply(valid_results, function(x) x$reject)
  stages <- sapply(valid_results, function(x) x$stage)
  decisions <- sapply(valid_results, function(x) x$decision)
  z1_vals <- sapply(valid_results, function(x) x$z1)
  z2_vals <- sapply(valid_results, function(x) x$z2)
  z_combo_vals <- sapply(valid_results, function(x) x$z_combo)
  w1_vals <- sapply(valid_results, function(x) x$w1)
  w2_vals <- sapply(valid_results, function(x) x$w2)

  empirical_rejection_rate <- mean(rejections)

  # Decision probabilities
  p_efficacy_emp <- mean(decisions == "efficacy")
  p_futility_emp <- mean(decisions == "futility")
  p_continue_emp <- mean(decisions == "continue")

  # Rejection by stage
  stage1_rejections <- sum(rejections & stages == 1, na.rm = TRUE)
  stage2_rejections <- sum(rejections & stages == 2, na.rm = TRUE)

  # Power breakdown
  power_s1 <- stage1_rejections / n_valid
  power_s2 <- stage2_rejections / n_valid

  # Theoretical values
  # Theoretical values
  # Theoretical values - extract scalars
  theoretical_power <- raw$power
  if (is.list(theoretical_power)) {
    theoretical_power <- theoretical_power$total %||% theoretical_power[[1]] %||% NA
  }
  if (length(theoretical_power) > 1) theoretical_power <- theoretical_power[1]

  probs <- raw$probabilities
  probs <- raw$probabilities

  # Extract scalar values safely
  p_eff <- probs$efficacy_stop
  if (length(p_eff) > 1) p_eff <- p_eff[1]

  p_fut <- probs$futility_stop
  if (length(p_fut) > 1) p_fut <- p_fut[1]

  p_cont <- probs$continue
  if (length(p_cont) > 1) p_cont <- p_cont[1]

  theoretical_power_s1 <- raw$power_stage1 %||% p_eff %||% NA
  theoretical_power_s2 <- raw$power_stage2 %||% NA

  if (length(theoretical_power_s2) == 1 && is.na(theoretical_power_s2) &&
      length(theoretical_power) == 1 && !is.na(theoretical_power) &&
      length(theoretical_power_s1) == 1 && !is.na(theoretical_power_s1)) {
    theoretical_power_s2 <- theoretical_power - theoretical_power_s1
  }

  # Update probs for later use
  probs$efficacy_stop <- p_eff
  probs$futility_stop <- p_fut
  probs$continue <- p_cont

  # Monte Carlo SE
  mc_se <- sqrt(empirical_rejection_rate * (1 - empirical_rejection_rate) / n_valid)

  if (verbose) {
    cat("\n=== Results ===\n")
    cat(sprintf("Valid simulations: %d / %d (%.1f%%)\n",
                n_valid, n_sims, 100 * n_valid / n_sims))

    cat(sprintf("\n--- %s ---\n", if (hypothesis == "H0") "Type I Error" else "Power"))
    cat(sprintf("Empirical: %.4f (%.2f%%)\n",
                empirical_rejection_rate, empirical_rejection_rate * 100))
    cat(sprintf("Theoretical: %.4f (%.2f%%)\n",
                if (hypothesis == "H0") 0.05 else theoretical_power,
                if (hypothesis == "H0") 5 else theoretical_power * 100))
    cat(sprintf("Monte Carlo SE: %.4f\n", mc_se))
    cat(sprintf("95%% CI: [%.4f, %.4f]\n",
                empirical_rejection_rate - 1.96 * mc_se,
                empirical_rejection_rate + 1.96 * mc_se))

    cat("\n--- Power by stage ---\n")
    cat(sprintf("Stage 1: Empirical=%.3f, Theoretical=%.3f\n",
                power_s1, theoretical_power_s1 %||% NA))
    cat(sprintf("Stage 2: Empirical=%.3f, Theoretical=%.3f\n",
                power_s2, theoretical_power_s2 %||% NA))

    cat("\n--- Decision probabilities ---\n")
    cat(sprintf("P(efficacy): Empirical=%.3f, Theoretical=%.3f\n",
                p_efficacy_emp, probs$efficacy_stop %||% NA))
    cat(sprintf("P(futility): Empirical=%.3f, Theoretical=%.3f\n",
                p_futility_emp, probs$futility_stop %||% NA))
    cat(sprintf("P(continue): Empirical=%.3f, Theoretical=%.3f\n",
                p_continue_emp, probs$continue %||% NA))

    cat("\n--- Z-statistic summaries ---\n")
    w1_ref <- stage2_weights[[names(stage2_weights)[1]]]$w1
    cat(sprintf("\nEffective stage 1 efficacy boundary: |z1| > %.3f (since w1=%.3f)\n",
                1.96/w1_ref, w1_ref))
    cat(sprintf("z1: mean=%.3f, sd=%.3f (expected under %s: mean=%.3f, sd=1)\n",
                mean(z1_vals, na.rm = TRUE), sd(z1_vals, na.rm = TRUE),
                hypothesis, if (hypothesis == "H0") 0 else raw$params$mu1 %||% NA))

    z2_continue <- z2_vals[!is.na(z2_vals)]
    if (length(z2_continue) > 0) {
      cat(sprintf("z2|continue: mean=%.3f, sd=%.3f\n",
                  mean(z2_continue), sd(z2_continue)))
    }

    z_combo_continue <- z_combo_vals[!is.na(z_combo_vals)]
    if (length(z_combo_continue) > 0) {
      cat(sprintf("z_combo|continue: mean=%.3f, sd=%.3f (expected sd≈1)\n",
                  mean(z_combo_continue), sd(z_combo_continue)))
    }

    # Add after Z-statistic summaries
    cat("\n--- Degrees of freedom (Satterthwaite) ---\n")
    df1_vals <- sapply(valid_results, function(x) x$df1 %||% NA)
    df2_vals <- sapply(valid_results, function(x) x$df2 %||% NA)
    cat(sprintf("Stage 1 df: mean=%.1f, range=[%.1f, %.1f]\n",
                mean(df1_vals, na.rm = TRUE),
                min(df1_vals, na.rm = TRUE),
                max(df1_vals, na.rm = TRUE)))
    df2_valid <- df2_vals[!is.na(df2_vals)]
    if (length(df2_valid) > 0) {
      cat(sprintf("Stage 2 df: mean=%.1f, range=[%.1f, %.1f]\n",
                  mean(df2_valid),
                  min(df2_valid),
                  max(df2_valid)))
    }

    cat("\n--- Raw t-statistics (before z-transform) ---\n")
    t1_vals <- sapply(valid_results, function(x) x$t1 %||% NA)
    t2_vals <- sapply(valid_results, function(x) x$t2 %||% NA)
    cat(sprintf("t1: mean=%.3f, sd=%.3f\n",
                mean(t1_vals, na.rm = TRUE), sd(t1_vals, na.rm = TRUE)))
    t2_valid <- t2_vals[!is.na(t2_vals)]
    if (length(t2_valid) > 0) {
      cat(sprintf("t2|continue: mean=%.3f, sd=%.3f\n",
                  mean(t2_valid), sd(t2_valid)))
    }

    # Check weights
    w1_continue <- w1_vals[!is.na(w1_vals)]
    w2_continue <- w2_vals[!is.na(w2_vals)]
    if (length(w1_continue) > 0) {
      cat(sprintf("\nWeights: w1 range=[%.3f, %.3f], w2 range=[%.3f, %.3f]\n",
                  min(w1_continue), max(w1_continue),
                  min(w2_continue), max(w2_continue)))
      cat(sprintf("w1^2 + w2^2 range: [%.3f, %.3f] (should be ≈1)\n",
                  min(w1_continue^2 + w2_continue^2),
                  max(w1_continue^2 + w2_continue^2)))
      z1_weighted_vals <- sapply(valid_results, function(x) {
        if (!is.null(x$w1) && !is.na(x$w1)) x$w1 * x$z1 else x$z1_weighted %||% NA
      })

      cat(sprintf("\nWeighted z1 (w1*z1): mean=%.3f, sd=%.3f\n",
                  mean(z1_weighted_vals, na.rm = TRUE),
                  sd(z1_weighted_vals, na.rm = TRUE)))
    }
  }

  list(
    empirical = list(
      rejection_rate = empirical_rejection_rate,
      mc_se = mc_se,
      power_s1 = power_s1,
      power_s2 = power_s2,
      p_efficacy = p_efficacy_emp,
      p_futility = p_futility_emp,
      p_continue = p_continue_emp,
      stage1_rejections = stage1_rejections,
      stage2_rejections = stage2_rejections,
      n_valid = n_valid
    ),
    theoretical = list(
      target = if (hypothesis == "H0") 0.05 else theoretical_power,
      power_s1 = theoretical_power_s1,
      power_s2 = theoretical_power_s2,
      p_efficacy = probs$efficacy_stop,
      p_futility = probs$futility_stop,
      p_continue = probs$continue,
      mu1 = raw$params$mu1
    ),
    z1 = z1_vals,
    z2 = z2_vals,
    z_combo = z_combo_vals,
    w1 = w1_vals,
    w2 = w2_vals,
    hypothesis = hypothesis,
    n_sims = n_sims,
    n_valid = n_valid,
    sim_results = valid_results
  )
}

#' Plot validation diagnostics
#' @param sim_output Output from simulate_adaptive_trial
plot_simulation_diagnostics <- function(sim_output) {
  library(ggplot2)

  hypothesis <- sim_output$hypothesis

  # Z1 distribution
  df_z1 <- data.frame(z = sim_output$z1[!is.na(sim_output$z1)])

  mu_theory <- if (hypothesis == "H0") 0 else sim_output$theoretical$mu1 %||% 0

  p1 <- ggplot(df_z1, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "#2a9d8f", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = mu_theory, sd = 1),
                  colour = "#e76f51", linewidth = 1) +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey40") +
    labs(x = expression(z[1]), y = "Density",
         title = expression("Stage 1: "*z[1]*" distribution"),
         subtitle = sprintf("Empirical: mean=%.2f, sd=%.2f | Theoretical: N(%.2f, 1)",
                            mean(df_z1$z), sd(df_z1$z), mu_theory)) +
    theme_minimal()

  # Z2 distribution (conditional on continuing)
  z2_valid <- sim_output$z2[!is.na(sim_output$z2)]

  if (length(z2_valid) > 10) {
    df_z2 <- data.frame(z = z2_valid)

    p2 <- ggplot(df_z2, aes(x = z)) +
      geom_histogram(aes(y = after_stat(density)), bins = 50,
                     fill = "#e9c46a", alpha = 0.7) +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                    colour = "#264653", linewidth = 1, linetype = "dashed") +
      labs(x = expression(z[2]), y = "Density",
           title = expression("Stage 2: "*z[2]*" | continue"),
           subtitle = sprintf("mean=%.2f, sd=%.2f (dashed: N(0,1) reference)",
                              mean(z2_valid), sd(z2_valid))) +
      theme_minimal()
  } else {
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient stage 2 data") +
      theme_void()
  }

  # Combination statistic
  z_combo_valid <- sim_output$z_combo[!is.na(sim_output$z_combo)]

  if (length(z_combo_valid) > 10) {
    df_combo <- data.frame(z = z_combo_valid)

    p3 <- ggplot(df_combo, aes(x = z)) +
      geom_histogram(aes(y = after_stat(density)), bins = 50,
                     fill = "#f4a261", alpha = 0.7) +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                    colour = "#264653", linewidth = 1, linetype = "dashed") +
      geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "#e76f51") +
      labs(x = expression(z[combo]), y = "Density",
           title = expression("Combination: "*z[combo]*" = "*w[1]*z[1]*" + "*w[2]*z[2]),
           subtitle = sprintf("mean=%.2f, sd=%.2f (expected: sd≈1)",
                              mean(z_combo_valid), sd(z_combo_valid))) +
      theme_minimal()
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient combination data") +
      theme_void()
  }

  # Decision comparison bar chart
  df_decisions <- data.frame(
    Decision = factor(rep(c("Efficacy", "Continue", "Futility"), 2),
                      levels = c("Efficacy", "Continue", "Futility")),
    Probability = c(
      sim_output$empirical$p_efficacy,
      sim_output$empirical$p_continue,
      sim_output$empirical$p_futility,
      sim_output$theoretical$p_efficacy %||% NA,
      sim_output$theoretical$p_continue %||% NA,
      sim_output$theoretical$p_futility %||% NA
    ),
    Source = factor(rep(c("Empirical", "Theoretical"), each = 3),
                    levels = c("Empirical", "Theoretical"))
  )

  p4 <- ggplot(df_decisions, aes(x = Decision, y = Probability, fill = Source)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Empirical" = "#2a9d8f", "Theoretical" = "#264653")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    labs(title = "Decision probabilities",
         subtitle = sprintf("N = %d simulations", sim_output$n_valid)) +
    theme_minimal() +
    theme(legend.position = "bottom")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    (p1 + p2) / (p3 + p4)
  } else {
    list(z1 = p1, z2 = p2, z_combo = p3, decisions = p4)
  }
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

#' Plot validation diagnostics
#' @param sim_output Output from simulate_adaptive_trial
plot_simulation_diagnostics <- function(sim_output) {

  hypothesis <- sim_output$hypothesis
  mu_theory <- if (hypothesis == "H0") 0 else sim_output$theoretical$target

  # Z1 distribution
  df_z1 <- data.frame(z = sim_output$z1[!is.na(sim_output$z1)])

  p1 <- ggplot(df_z1, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = "#2a9d8f", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  colour = "#264653", linewidth = 1) +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "#e76f51") +
    labs(x = expression(z[1]), y = "Density",
         title = expression("Stage 1: "*z[1]*" distribution"),
         subtitle = sprintf("Mean=%.2f, SD=%.2f (expected: SD=1)",
                            mean(df_z1$z), sd(df_z1$z))) +
    theme_minimal()

  # Z2 distribution (conditional on continuing)
  z2_valid <- sim_output$z2[!is.na(sim_output$z2)]

  if (length(z2_valid) > 10) {
    df_z2 <- data.frame(z = z2_valid)

    p2 <- ggplot(df_z2, aes(x = z)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                     fill = "#e9c46a", alpha = 0.7) +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                    colour = "#264653", linewidth = 1) +
      labs(x = expression(z[2]), y = "Density",
           title = expression("Stage 2: "*z[2]*" | continue"),
           subtitle = sprintf("Mean=%.2f, SD=%.2f", mean(z2_valid), sd(z2_valid))) +
      theme_minimal()
  } else {
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Not enough stage 2 data") +
      theme_void()
  }

  # Combination statistic
  z_combo_valid <- sim_output$z_combo[!is.na(sim_output$z_combo)]

  if (length(z_combo_valid) > 10) {
    df_combo <- data.frame(z = z_combo_valid)

    p3 <- ggplot(df_combo, aes(x = z)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                     fill = "#f4a261", alpha = 0.7) +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                    colour = "#264653", linewidth = 1) +
      geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "#e76f51") +
      labs(x = expression(z[combo]), y = "Density",
           title = expression("Combination: "*z[combo]*" | continue"),
           subtitle = sprintf("Mean=%.2f, SD=%.2f (expected: SD=1)",
                              mean(z_combo_valid), sd(z_combo_valid))) +
      theme_minimal()
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Not enough combination data") +
      theme_void()
  }

  # Decision comparison
  df_decisions <- data.frame(
    Decision = rep(c("Efficacy", "Futility", "Continue"), 2),
    Probability = c(
      sim_output$empirical$p_efficacy,
      sim_output$empirical$p_futility,
      sim_output$empirical$p_continue,
      sim_output$theoretical$p_efficacy,
      sim_output$theoretical$p_futility,
      sim_output$theoretical$p_continue
    ),
    Source = rep(c("Empirical", "Theoretical"), each = 3)
  )

  p4 <- ggplot(df_decisions, aes(x = Decision, y = Probability, fill = Source)) +
    geom_col(position = "dodge", width = 0.6) +
    scale_fill_manual(values = c("Empirical" = "#2a9d8f", "Theoretical" = "#264653")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Decision probabilities", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "bottom")

  (p1 + p2) / (p3 + p4)
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
