# =============================================================================
# Simplified Analysis Interface
# =============================================================================
# Provides a clean, user-friendly interface for running adaptive trial analyses.
# Users work with design objects and get results objects with built-in methods.

# =============================================================================
# Main Analysis Function
# =============================================================================

#' Run adaptive trial analysis
#'
#' This is the main entry point for analysing adaptive CRT designs. It handles
#' all the complexity of lambda calibration, design space exploration, and
#' result compilation.
#'
#' @param design A crt_design object created by parallel_crt(), crossover_crt(), etc.
#' @param target_power Target overall power (default 0.8)
#' @param explore If TRUE and stage1_grid has multiple rows, explore all designs.
#'        If FALSE, only analyse the first stage 1 design.
#' @param verbose Print progress messages
#' @param ... Additional arguments passed to find_lambda_for_power
#'
#' @return An adaptive_results object with methods for plotting and summarising
#'
#' @examples
#' design <- parallel_crt(
#'   icc = 0.05, delta = 0.25,
#'   k1 = c(8, 10, 12), m1 = seq(10, 50, by = 10),
#'   k2 = 0:4, rho = 30
#' )
#' results <- adaptive_analysis(design, target_power = 0.8)
#' summary(results)
#' plot(results, type = "EN")
#'
#' @export
adaptive_analysis <- function(design, target_power = 0.8,
                              explore = TRUE, verbose = TRUE, ...) {
  UseMethod("adaptive_analysis")
}

#' @export
adaptive_analysis.crt_design <- function(design, target_power = 0.8,
                                         explore = TRUE, verbose = TRUE, ...) {

  spec <- design$spec
  rho <- design$rho

  # Generate functions from specification
  cost_fn <- generate_cost_fn(spec, "stage2")
  sample_size_fn <- generate_sample_size_fn(spec)
  model_fn_factory <- generate_model_factory(spec)

  # Decide whether to explore or single analysis
  do_explore <- explore && nrow(design$stage1_grid) > 1

  if (do_explore) {
    if (verbose) {
      cat(sprintf("Exploring %d stage 1 designs for %.0f%% power...\n\n",
                  nrow(design$stage1_grid), target_power * 100))
    }

    raw_results <- explore_design_space(
      stage1_grid = design$stage1_grid,
      target_power = target_power,
      design_grid_fn = design$stage2_grid_fn,
      model_fn_factory = model_fn_factory,
      fixed_params_base = design$fixed_params,
      cost_fn = cost_fn,
      cost_params = list(rho = rho),
      resource_vars = spec$stage2_params,
      sample_size_fn = sample_size_fn,
      verbose = verbose,
      ...
    )

    # Add sample size columns if missing (fallback)
    if (!"E_N" %in% names(raw_results$summary)) {
      raw_results <- add_sample_size_metrics(raw_results, sample_size_fn)
    }

  } else {
    # Single design analysis
    s1 <- design$stage1_grid[1, , drop = FALSE]
    fixed_params <- c(design$fixed_params, as.list(s1))

    if (verbose) {
      param_str <- paste(names(s1), "=", s1, collapse = ", ")
      cat(sprintf("Analysing single design: %s\n", param_str))
    }

    raw_results <- find_lambda_for_power(
      target_power = target_power,
      design_grid = design$stage2_grid_fn(s1),
      model_fn = model_fn_factory(),
      fixed_params = fixed_params,
      cost_fn = cost_fn,
      cost_params = list(rho = rho),
      resource_vars = spec$stage2_params,
      verbose = verbose,
      ...
    )
  }

  structure(
    list(
      raw = raw_results,
      design = design,
      target_power = target_power,
      explored = do_explore,
      sample_size_fn = sample_size_fn
    ),
    class = "adaptive_results"
  )
}

#' Add sample size metrics to exploration results (fallback)
#' @keywords internal
add_sample_size_metrics <- function(exploration, sample_size_fn) {

  for (i in seq_len(nrow(exploration$summary))) {
    if (!is.null(exploration$results[[i]]$results)) {
      ss <- tryCatch({
        compute_sample_size(
          exploration$results[[i]]$results,
          sample_size_fn,
          hypothesis = "H1"
        )
      }, error = function(e) NULL)

      if (!is.null(ss)) {
        for (j in seq_len(nrow(ss$summary))) {
          col_name <- ss$summary$metric[j]
          if (i == 1 && !col_name %in% names(exploration$summary)) {
            exploration$summary[[col_name]] <- NA
          }
          exploration$summary[[col_name]][i] <- ss$summary$value[j]
        }
      }
    }
  }

  exploration
}

# =============================================================================
# Methods for adaptive_results
# =============================================================================

#' @export
print.adaptive_results <- function(x, ...) {
  cat("Adaptive CRT Analysis Results\n")
  cat("=============================\n\n")

  cat(sprintf("Design type: %s\n", x$design$spec$design_type))
  cat(sprintf("Target power: %.0f%%\n", x$target_power * 100))

  if (x$explored) {
    cat(sprintf("Stage 1 designs explored: %d\n", nrow(x$raw$summary)))
    cat(sprintf("Converged: %d / %d\n",
                sum(x$raw$summary$converged, na.rm = TRUE),
                nrow(x$raw$summary)))
  } else {
    cat(sprintf("Converged: %s\n", x$raw$converged))
    cat(sprintf("Achieved power: %.3f\n", x$raw$power))
    cat(sprintf("Lambda: %.4e\n", x$raw$lambda))
    if (!is.null(x$raw$results$params$efficacy_boundary)) {
      cat(sprintf("Efficacy boundary: |z1| > %.3f\n",
                  x$raw$results$params$efficacy_boundary))
    }
  }

  cat("\nUse summary() for detailed results, plot() for visualisation.\n")
  invisible(x)
}

#' @export
summary.adaptive_results <- function(object, ...) {
  if (object$explored) {
    # Return exploration summary
    df <- object$raw$summary

    # Select key columns
    key_cols <- c(
      names(object$design$spec$stage1_params),
      "converged", "power", "power_stage1", "lambda",
      "E_N", "N_max", "E_cost", "max_cost",
      "prob_efficacy", "prob_futility", "prob_continue"
    )
    key_cols <- intersect(key_cols, names(df))

    cat("Design Space Exploration Summary\n")
    cat("================================\n\n")

    print(df[, key_cols], row.names = FALSE)
    invisible(df)

  } else {
    # Single design summary
    cat("Single Design Analysis Summary\n")
    cat("==============================\n\n")

    cat(sprintf("Converged: %s\n", object$raw$converged))
    cat(sprintf("Target power: %.3f\n", object$target_power))
    cat(sprintf("Achieved power: %.3f\n", object$raw$power))
    cat(sprintf("Lambda: %.4e\n\n", object$raw$lambda))

    if (!is.null(object$raw$results)) {
      cat("Power breakdown:\n")
      cat(sprintf("  Stage 1: %.3f\n", object$raw$results$power$stage1))
      cat(sprintf("  Stage 2: %.3f\n", object$raw$results$power$stage2))

      cat("\nDecision probabilities:\n")
      cat(sprintf("  Efficacy stop: %.3f\n", object$raw$results$probabilities$efficacy_stop))
      cat(sprintf("  Futility stop: %.3f\n", object$raw$results$probabilities$futility_stop))
      cat(sprintf("  Continue: %.3f\n", object$raw$results$probabilities$continue))

      # Compute sample sizes
      ss <- compute_sample_size(object$raw$results, object$sample_size_fn, "H1")
      cat("\nSample size summary:\n")
      print(ss$summary, row.names = FALSE)
    }

    invisible(object$raw)
  }
}

#' @export
plot.adaptive_results <- function(x, type = c("EN", "Nmax", "power", "probabilities",
                                              "decision", "pareto", "all"),
                                  objectives = NULL,
                                  design_index = 1,
                                  ...) {
  type <- match.arg(type)

  if (x$explored) {
    # Exploration plots
    if (type == "decision") {
      # Show decision rules for specified design
      idx <- design_index
      if (!is.null(x$raw$results[[idx]]$results)) {
        plot_decision_rules(x, design_index = idx, ...)
      } else {
        message("No converged results available for decision plot")
      }
    } else if (type == "pareto") {
      # Pareto frontier - use provided objectives or default
      if (is.null(objectives)) {
        objectives <- list(E_N = "min", N_max = "min")
      }
      pareto <- find_constrained_pareto(x$raw$summary, objectives = objectives)
      plot_pareto(pareto, all_designs = x$raw$summary, objectives = objectives, ...)
    } else {
      # Standard exploration plots
      plot_design_space(x$raw, type = type, ...)
    }
  } else {
    # Single design plots
    if (type == "decision") {
      plot_decision_rules(x, design_index = 1, ...)
    } else {
      message("For single designs, use type = 'decision'")
    }
  }
}

# =============================================================================
# Pareto Frontier for Results
# =============================================================================

#' Find Pareto frontier from adaptive analysis results
#'
#' @param results An adaptive_results object
#' @param objectives Named list of objectives (default: E_N and N_max minimised)
#' @param constraints Optional constraints (see find_constrained_pareto)
#' @param ... Additional arguments passed to find_constrained_pareto
#'
#' @return A data frame of Pareto-optimal designs
#' @export
find_pareto <- function(results,
                        objectives = list(E_N = "min", N_max = "min"),
                        constraints = NULL, ...) {

  if (!inherits(results, "adaptive_results")) {
    stop("results must be an adaptive_results object")
  }

  if (!results$explored) {
    stop("Pareto analysis requires explored results (multiple stage 1 designs)")
  }

  find_constrained_pareto(results$raw$summary, objectives, constraints, ...)
}

# =============================================================================
# Fixed Design Comparison
# =============================================================================

#' Compare adaptive design to optimal fixed design
#'
#' @param results An adaptive_results object
#' @param fixed_grid Data frame of fixed design parameters to search over.
#'        If NULL, auto-generated based on design type.
#' @param fixed_model_fn Custom fixed model function. If NULL, uses design-type default.
#' @param fixed_cost_fn Custom fixed cost function. If NULL, uses design-type default.
#' @param verbose Print progress
#'
#' @return A comparison object with comparison metrics
#' @export
compare_to_fixed <- function(results,
                             fixed_grid = NULL,
                             fixed_model_fn = NULL,
                             fixed_cost_fn = NULL,
                             verbose = TRUE) {

  if (!inherits(results, "adaptive_results")) {
    stop("results must be an adaptive_results object")
  }

  design <- results$design
  spec <- design$spec
  rho <- design$rho

  # Get design-type-specific defaults if not provided
  fixed_defaults <- get_fixed_design_defaults(design, results$target_power)

  if (is.null(fixed_grid)) {
    fixed_grid <- fixed_defaults$grid
    if (is.null(fixed_grid)) {
      stop("Cannot auto-generate fixed design grid for design type '", spec$design_type,
           "'. Please provide fixed_grid.")
    }
    if (verbose) {
      cat("Using default fixed design grid for", spec$design_type, "design:\n")
      cat("  Parameters:", paste(names(fixed_grid), collapse = ", "), "\n")
      cat("  Grid size:", nrow(fixed_grid), "designs\n\n")
    }
  }

  if (is.null(fixed_model_fn)) {
    fixed_model_fn <- fixed_defaults$model_fn
    if (is.null(fixed_model_fn)) {
      stop("No default fixed model function for design type '", spec$design_type,
           "'. Please provide fixed_model_fn.")
    }
  }

  if (is.null(fixed_cost_fn)) {
    fixed_cost_fn <- fixed_defaults$cost_fn
    if (is.null(fixed_cost_fn)) {
      fixed_cost_fn <- function(model_output, rho) {
        model_output$total_cost %||% model_output$n_total %||% 0
      }
    }
  }

  if (verbose) cat("Finding optimal fixed design...\n")

  fixed_result <- find_optimal_fixed_design(
    target_power = results$target_power,
    design_grid = fixed_grid,
    model_fn_factory = function() fixed_model_fn,
    fixed_params_base = design$fixed_params,
    cost_fn = fixed_cost_fn,
    cost_params = list(rho = rho),
    verbose = verbose
  )

  # Check if we found a valid design
  if (is.null(fixed_result$optimal) || length(fixed_result$optimal$cost) == 0) {
    warning("No fixed design achieved target power. Try expanding fixed_grid.\n",
            "  Max power achieved: ",
            sprintf("%.1f%%", max(fixed_result$all$power, na.rm = TRUE) * 100),
            "\n  Consider increasing cluster count or sample size.")

    # Return comparison with NA values
    if (results$explored) {
      comparison_df <- results$raw$summary
      comparison_df$fixed_cost <- NA
      comparison_df$fixed_power <- NA
      comparison_df$E_cost_ratio <- NA
      comparison_df$max_cost_ratio <- NA
      comparison_df$E_savings_pct <- NA
      comparison_df$max_within_fixed <- NA
    } else {
      comparison_df <- data.frame(
        E_cost = NA,
        max_cost = NA,
        fixed_cost = NA,
        E_cost_ratio = NA,
        max_cost_ratio = NA,
        E_savings_pct = NA,
        max_within_fixed = NA
      )
    }

    return(structure(
      list(
        comparison = comparison_df,
        fixed = fixed_result,
        adaptive = results,
        stage1_grid = design$stage1_grid,
        rho = rho,
        valid = FALSE
      ),
      class = "adaptive_comparison"
    ))
  }

  if (verbose) {
    cat("\nOptimal fixed design:\n")
    print(fixed_result$optimal)
  }

  # Build comparison data
  if (results$explored) {
    df <- results$raw$summary
    fixed_cost <- fixed_result$optimal$cost

    # Ensure cost columns exist
    if (!"E_cost" %in% names(df)) {
      df <- add_sample_size_metrics_to_df(df, results)
    }

    if (!"E_cost" %in% names(df)) {
      warning("E_cost not found; using E_N as proxy")
      df$E_cost <- df$E_N
      df$max_cost <- df$N_max
    }

    df$fixed_cost <- fixed_cost
    df$fixed_power <- fixed_result$optimal$power
    df$E_cost_ratio <- df$E_cost / fixed_cost
    df$max_cost_ratio <- df$max_cost / fixed_cost
    df$E_savings_pct <- (1 - df$E_cost_ratio) * 100
    df$max_within_fixed <- df$max_cost <= fixed_cost

    comparison_df <- df
  } else {
    # Single design comparison
    ss <- compute_sample_size(results$raw$results, results$sample_size_fn, "H1")
    E_cost <- ss$summary$value[ss$summary$metric == "E_cost"]
    max_cost <- ss$summary$value[ss$summary$metric == "max_cost"]

    if (length(E_cost) == 0) E_cost <- ss$summary$value[ss$summary$metric == "E_N"]
    if (length(max_cost) == 0) max_cost <- ss$summary$value[ss$summary$metric == "N_max"]

    fixed_cost <- fixed_result$optimal$cost

    comparison_df <- data.frame(
      E_cost = E_cost,
      max_cost = max_cost,
      fixed_cost = fixed_cost,
      E_cost_ratio = E_cost / fixed_cost,
      max_cost_ratio = max_cost / fixed_cost,
      E_savings_pct = (1 - E_cost / fixed_cost) * 100,
      max_within_fixed = max_cost <= fixed_cost
    )
  }

  structure(
    list(
      comparison = comparison_df,
      fixed = fixed_result,
      adaptive = results,
      stage1_grid = design$stage1_grid,
      rho = rho,
      valid = TRUE
    ),
    class = "adaptive_comparison"
  )
}

#' Get default fixed design configuration for a design type
#' @keywords internal
get_fixed_design_defaults <- function(design, target_power) {

  spec <- design$spec
  rho <- design$rho
  stage1_grid <- design$stage1_grid
  fixed_params <- design$fixed_params

  if (spec$design_type == "parallel") {
    # Parallel: vary k and m - use wider range
    k_vals <- unique(stage1_grid$k1)
    m_vals <- unique(stage1_grid$m1)

    # Expand grid more generously to ensure we can find a design with target power
    k_max <- max(ceiling(max(k_vals) * 2.5), max(k_vals) + 20)
    m_max <- max(ceiling(max(m_vals) * 3), max(m_vals) + 100)

    grid <- expand.grid(
      k = seq(min(k_vals), k_max, by = 2),
      m = seq(min(m_vals), m_max, by = 10)
    )

    model_fn <- fixed_parallel_model_builder

    cost_fn <- function(model_output, rho) {
      n_cost <- model_output$n_total
      k_cost <- model_output$k * 2
      n_cost + rho * k_cost
    }

  } else if (spec$design_type == "crossover") {
    # Crossover: vary k and m
    k_vals <- unique(stage1_grid$k1)
    m_vals <- unique(stage1_grid$m1)

    k_max <- max(ceiling(max(k_vals) * 2.5), max(k_vals) + 20)
    m_max <- max(ceiling(max(m_vals) * 3), max(m_vals) + 100)

    grid <- expand.grid(
      k = seq(min(k_vals), k_max, by = 2),
      m = seq(min(m_vals), m_max, by = 10)
    )

    model_fn <- fixed_crossover_model_builder

    cost_fn <- function(model_output, rho) {
      n_cost <- model_output$n_total
      k_cost <- model_output$n_clusters
      n_cost + rho * k_cost
    }
    } else {
    # Unknown design type - return NULL, user must provide
    return(list(
      grid = NULL,
      model_fn = NULL,
      cost_fn = NULL
    ))
  }

  list(
    grid = grid,
    model_fn = model_fn,
    cost_fn = cost_fn
  )
}

#' Helper to add sample size metrics to summary df
#' @keywords internal
add_sample_size_metrics_to_df <- function(df, results) {
  for (i in seq_len(nrow(df))) {
    if (!is.null(results$raw$results[[i]]$results)) {
      ss <- tryCatch({
        compute_sample_size(
          results$raw$results[[i]]$results,
          results$sample_size_fn,
          hypothesis = "H1"
        )
      }, error = function(e) NULL)

      if (!is.null(ss)) {
        for (j in seq_len(nrow(ss$summary))) {
          col_name <- ss$summary$metric[j]
          if (!col_name %in% names(df)) df[[col_name]] <- NA
          df[[col_name]][i] <- ss$summary$value[j]
        }
      }
    }
  }
  df
}

#' @export
print.adaptive_comparison <- function(x, ...) {
  cat("Adaptive vs Fixed Design Comparison\n")
  cat("===================================\n\n")

  if (!isTRUE(x$valid)) {
    cat("WARNING: No fixed design achieved target power.\n")
    if (!is.null(x$fixed$all)) {
      cat(sprintf("  Max power in grid: %.1f%%\n",
                  max(x$fixed$all$power, na.rm = TRUE) * 100))
    }
    cat("  Expand fixed_grid and re-run.\n\n")
    return(invisible(x))
  }

  # Print fixed design info
  fixed_opt <- x$fixed$optimal
  fixed_params <- setdiff(names(fixed_opt), c("power", "cost"))
  param_str <- paste(sapply(fixed_params, function(p) sprintf("%s = %g", p, fixed_opt[[p]])),
                     collapse = ", ")

  cat(sprintf("Fixed design: %s\n", param_str))
  cat(sprintf("Fixed power: %.3f\n", fixed_opt$power))
  cat(sprintf("Fixed cost: %.1f\n\n", fixed_opt$cost))

  if ("comparison" %in% names(x) && is.data.frame(x$comparison)) {
    if (nrow(x$comparison) > 1) {
      cat("Summary of adaptive designs:\n")
      cat(sprintf("  Designs evaluated: %d\n", nrow(x$comparison)))
      cat(sprintf("  E[cost] ratio range: %.2f - %.2f\n",
                  min(x$comparison$E_cost_ratio, na.rm = TRUE),
                  max(x$comparison$E_cost_ratio, na.rm = TRUE)))
      cat(sprintf("  Designs with max cost <= fixed: %d\n",
                  sum(x$comparison$max_within_fixed, na.rm = TRUE)))
    } else {
      cat(sprintf("E[cost] ratio: %.2f\n", x$comparison$E_cost_ratio))
      cat(sprintf("Max cost ratio: %.2f\n", x$comparison$max_cost_ratio))
      cat(sprintf("Expected savings: %.1f%%\n", x$comparison$E_savings_pct))
    }
  }

  invisible(x)
}

#' @export
plot.adaptive_comparison <- function(x, type = c("cost_ratio", "savings", "frontier", "all"),
                                     ...) {
  type <- match.arg(type)

  df <- x$comparison
  df <- df[df$converged,]

  if (!is.data.frame(df) || nrow(df) == 0) {
    message("No comparison data available for plotting")
    return(invisible(NULL))
  }

  if (nrow(df) == 1) {
    # Single design - show a simple bar comparison


    plot_data <- data.frame(
      metric = c("E[cost]", "Max cost", "Fixed cost"),
      value = c(df$E_cost, df$max_cost, df$fixed_cost),
      type = c("Adaptive", "Adaptive", "Fixed")
    )

    p <- ggplot(plot_data, aes(x = metric, y = value, fill = type)) +
      geom_col(position = "dodge", width = 0.6) +
      scale_fill_manual(values = c("Adaptive" = "#2a9d8f", "Fixed" = "#e76f51")) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = NULL, y = "Cost", fill = NULL,
           title = "Adaptive vs Fixed Design Cost",
           subtitle = sprintf("Expected savings: %.1f%%", df$E_savings_pct)) +
      theme_publication() +
      theme(legend.position = "bottom")

    return(p)
  }

  # Multiple designs - use custom plot function
  plot_comparison_multi(x, type = type, ...)
}

#' Plot comparison for multiple adaptive designs
#' @keywords internal
plot_comparison_multi <- function(x, type = c("cost_ratio", "savings", "frontier", "all")) {

  type <- match.arg(type)

  df <- x$comparison
  fixed_cost <- x$fixed$optimal$cost
  stage1_grid <- x$stage1_grid

  # Split by convergence instead of filtering
  df_conv <- df[df$converged, ]
  df_fail <- df[!df$converged, ]

  # Determine grouping variable
  stage1_vars <- names(stage1_grid)
  if (length(stage1_vars) == 0) {
    stage1_vars <- intersect(c("k1", "m1", "K", "t1"), names(df))
  }

  if (length(stage1_vars) == 0) {
    message("Cannot determine stage 1 variables for grouping")
    return(invisible(NULL))
  }

  group_var <- stage1_vars[1]

  if (!group_var %in% names(df)) {
    message(sprintf("Group variable '%s' not found in comparison data", group_var))
    return(invisible(NULL))
  }

  # Factor the group var in both subsets (same levels)
  all_levels <- sort(unique(df[[group_var]]))
  df_conv[[group_var]] <- factor(df_conv[[group_var]], levels = all_levels)
  df_fail[[group_var]] <- factor(df_fail[[group_var]], levels = all_levels)

  n_groups <- length(all_levels)
  palette_groups <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:min(n_groups, 5)]
  if (n_groups > 5) palette_groups <- scales::hue_pal()(n_groups)

  has_failures <- nrow(df_fail) > 0

  # Determine x variable
  x_var <- if ("power_stage1" %in% names(df)) {
    "power_stage1"
  } else if (length(stage1_vars) > 1) {
    stage1_vars[2]
  } else {
    group_var
  }

  if (!x_var %in% names(df)) x_var <- group_var

  # Helper to add non-converged point layer
  add_fail_layer <- function(p, y_var) {
    if (has_failures) {
      p <- p +
        geom_point(data = df_fail, aes(x = .data[[x_var]], y = .data[[y_var]]),
                   colour = "grey45", shape = 4, size = 2.5, stroke = 0.8,
                   inherit.aes = FALSE)
    }
    p
  }

  # Cost ratio plot
  p_ratio <- ggplot(df_conv, aes(x = .data[[x_var]])) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(aes(y = E_cost_ratio, colour = .data[[group_var]]), linewidth = 0.9) +
    geom_line(aes(y = max_cost_ratio, colour = .data[[group_var]]),
              linetype = "dotted", linewidth = 0.7) +
    geom_point(aes(y = E_cost_ratio, colour = .data[[group_var]]), size = 2) +
    scale_colour_manual(values = palette_groups) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = if (x_var == "power_stage1") "Stage 1 power" else gsub("_", " ", x_var),
      y = "Cost ratio (adaptive / fixed)",
      colour = group_var,
      title = "Cost comparison: adaptive vs fixed",
      subtitle = "Solid = E[cost], Dotted = max cost"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  p_ratio <- add_fail_layer(p_ratio, "E_cost_ratio")

  if (x_var == "power_stage1") {
    p_ratio <- p_ratio + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  }

  # Savings plot
  p_savings <- ggplot(df_conv, aes(x = .data[[x_var]], y = E_savings_pct)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(aes(colour = .data[[group_var]]), linewidth = 0.9) +
    geom_point(aes(colour = .data[[group_var]],
                   shape = factor(max_within_fixed, levels = c(TRUE, FALSE))), size = 3) +
    scale_colour_manual(values = palette_groups) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                       labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       name = "Max <= fixed") +
    labs(
      x = if (x_var == "power_stage1") "Stage 1 power" else gsub("_", " ", x_var),
      y = "Expected savings (%)",
      colour = group_var,
      title = "Expected cost savings vs fixed design"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  p_savings <- add_fail_layer(p_savings, "E_savings_pct")

  if (x_var == "power_stage1") {
    p_savings <- p_savings + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  }

  # Frontier plot
  p_frontier <- ggplot(df_conv, aes(x = E_cost, y = max_cost)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = fixed_cost, linetype = "dashed", colour = "#e76f51") +
    geom_hline(yintercept = fixed_cost, linetype = "dashed", colour = "#e76f51") +
    geom_point(aes(colour = .data[[group_var]],
                   shape = factor(max_within_fixed, levels = c(TRUE, FALSE))), size = 3) +
    scale_colour_manual(values = palette_groups) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                       labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       name = "Max <= fixed") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x = "E[cost]",
      y = "Max cost",
      colour = group_var,
      title = "Cost trade-off frontier",
      subtitle = "Dashed lines = fixed design cost"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  if (has_failures) {
    p_frontier <- p_frontier +
      geom_point(data = df_fail, aes(x = E_cost, y = max_cost),
                 colour = "grey45", shape = 4, size = 2.5, stroke = 0.8,
                 inherit.aes = FALSE)
  }

  if (type == "cost_ratio") return(p_ratio)
  if (type == "savings") return(p_savings)
  if (type == "frontier") return(p_frontier)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    (p_ratio + p_savings) / p_frontier +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    list(cost_ratio = p_ratio, savings = p_savings, frontier = p_frontier)
  }
}

# =============================================================================
# Sample Size Extraction
# =============================================================================

#' Get sample size distribution from results
#'
#' @param results An adaptive_results object
#' @param design_index For exploration results, which design to use (default: 1)
#' @param hypothesis "H1" or "H0"
#'
#' @return A sample size summary object
#' @export
get_sample_size <- function(results, design_index = 1, hypothesis = "H1") {

  if (!inherits(results, "adaptive_results")) {
    stop("results must be an adaptive_results object")
  }

  if (results$explored) {
    res <- results$raw$results[[design_index]]$results
  } else {
    res <- results$raw$results
  }

  if (is.null(res)) {
    stop("No results available for sample size computation")
  }

  compute_sample_size(res, results$sample_size_fn, hypothesis)
}

# =============================================================================
# Decision Rules Extraction
# =============================================================================

#' Get decision rules from results
#'
#' @param results An adaptive_results object
#' @param design_index For exploration results, which design to use
#' @param z1_values Optional specific z1 values to evaluate
#'
#' @return Data frame of decision rules
#' @export
get_decision_rules <- function(results, design_index = 1, z1_values = NULL) {

  if (!inherits(results, "adaptive_results")) {
    stop("results must be an adaptive_results object")
  }

  if (results$explored) {
    res <- results$raw$results[[design_index]]$results
  } else {
    res <- results$raw$results
  }

  if (is.null(res)) {
    stop("No results available")
  }

  opt <- res$quadrature$optimal_designs

  if (!is.null(z1_values)) {
    # Interpolate to specific z1 values
    # For now, just return nearest points
    idx <- sapply(z1_values, function(z) which.min(abs(opt$z1 - z)))
    opt <- opt[idx, ]
  }

  opt
}

#' Perform interim analysis with updated parameter estimates
#'
#' Uses pre-planned combination test weights but re-evaluates stage 2
#' designs using interim estimates of auxiliary parameters.
#'
#' @param planned An adaptive_results object from adaptive_analysis()
#' @param z1_obs Observed stage 1 test statistic
#' @param theta_hat Named list of updated parameter estimates (e.g. list(icc = 0.03, sigma2 = 1.2))
#' @param delta Treatment effect for conditional power (NULL = use planned value)
#' @param verbose Print details
#'
#' @return An interim_result object with recommendation and conditional power
interim_analysis <- function(planned, z1_obs, theta_hat = list(),
                             delta = NULL, recalibrate = TRUE,
                             target_power = 0.8, verbose = TRUE) {

  if (!inherits(planned, "adaptive_results")) {
    stop("planned must be an adaptive_results object from adaptive_analysis()")
  }

  design <- planned$design
  spec <- design$spec
  rho <- design$rho

  if (planned$explored) {
    stop("planned must be a single-design result. ",
         "Use adaptive_analysis() with a single stage 1 design, or extract ",
         "one result from the explored set.")
  }

  r <- planned$raw
  w1_ref <- r$results$params$w1_ref
  efficacy_boundary <- r$results$params$efficacy_boundary
  t_crit <- r$results$params$t_crit %||% 1.96
  b1_planned <- r$results$params$b1
  method <- r$method %||% "lambda"

  if (is.null(delta)) delta <- b1_planned

  if (verbose) {
    cat("=== Interim Analysis ===\n")
    cat(sprintf("Observed z1 = %.3f\n", z1_obs))
    cat(sprintf("Pre-planned weights: w1 = %.4f\n", w1_ref))
    cat(sprintf("Efficacy boundary: |z1| > %.3f\n", efficacy_boundary))
    if (length(theta_hat) > 0) {
      cat("Updated parameters:\n")
      for (nm in names(theta_hat)) {
        cat(sprintf("  %s = %s\n", nm, format(theta_hat[[nm]], digits = 4)))
      }
    }
  }

  # --- Step 1: Check efficacy stopping ---
  if (abs(z1_obs) > efficacy_boundary) {
    if (verbose) cat("\n** RECOMMENDATION: Stop for efficacy **\n")
    return(structure(list(
      decision = "efficacy_stop",
      z1 = z1_obs,
      w1_ref = w1_ref,
      weighted_statistic = w1_ref * z1_obs,
      critical_value = t_crit,
      stage2_design = NULL,
      conditional_power = NA,
      message = sprintf("z1 = %.3f exceeds efficacy boundary %.3f",
                        abs(z1_obs), efficacy_boundary)
    ), class = "interim_result"))
  }

  # --- Step 2: Rebuild models with updated parameters ---
  if (verbose) cat("\nRecomputing stage 2 designs with updated parameters...\n")

  fixed_params <- design$fixed_params
  s1 <- design$stage1_grid[1, , drop = FALSE]
  for (nm in names(s1)) fixed_params[[nm]] <- s1[[nm]]
  for (nm in names(theta_hat)) fixed_params[[nm]] <- theta_hat[[nm]]

  design_grid <- design$stage2_grid_fn(s1)
  model_fn <- generate_model_factory(spec)()
  cost_fn <- generate_cost_fn(spec, "stage2")
  resource_vars <- spec$stage2_params

  # --- Step 3: Recalibrate or use planned lambda ---
  if (recalibrate && length(theta_hat) > 0) {
    if (verbose) cat("Recalibrating lambda with updated parameters (weights locked)...\n")

    recal <- find_lambda_for_power(
      target_power = target_power,
      design_grid = design_grid,
      model_fn = model_fn,
      fixed_params = fixed_params,
      cost_fn = cost_fn,
      cost_params = list(rho = rho),
      resource_vars = resource_vars,
      method = method,
      verbose = verbose,
      w1_override = w1_ref,
      efficacy_override = efficacy_boundary,
      t_crit_override = t_crit
    )

    lambda_use <- recal$lambda
    cost_cap_use <- recal$cost_cap

    if (verbose) {
      cat(sprintf("Planned lambda: %.4e\n", r$lambda %||% NA))
      cat(sprintf("Recalibrated lambda: %.4e\n", lambda_use %||% NA))
    }
  } else {
    lambda_use <- r$lambda
    cost_cap_use <- r$cost_cap
  }

  # --- Step 4: Compute CP for each design ---
  model_list <- precompute_models(design_grid, model_fn, fixed_params,
                                  parallel = FALSE, n_cores = NULL)

  n_designs <- length(model_list)
  w2_ref <- sqrt(1 - w1_ref^2)

  cp_vec <- numeric(n_designs)
  cost_vec <- numeric(n_designs)
  I2_eff_vec <- numeric(n_designs)

  for (j in seq_len(n_designs)) {
    mod <- model_list[[j]]
    I2_eff_j <- mod$I2_eff

    ncp <- delta * sqrt(I2_eff_j)
    upper_thresh <- (t_crit - w1_ref * z1_obs) / w2_ref
    lower_thresh <- (-t_crit - w1_ref * z1_obs) / w2_ref

    cp_vec[j] <- pnorm(ncp - upper_thresh) + pnorm(lower_thresh - ncp)
    I2_eff_vec[j] <- I2_eff_j
    cost_vec[j] <- do.call(cost_fn, c(list(mod), list(rho = rho)))
  }

  # --- Step 5: Apply decision criterion ---
  if (method == "cost_cap") {
    feasible <- cost_vec <= cost_cap_use

    if (verbose) cat(sprintf("Budget constraint: cost <= %.1f\n", cost_cap_use))

    if (!any(feasible)) {
      if (verbose) cat("\n** RECOMMENDATION: Stop for futility (no feasible designs) **\n")
      return(structure(list(
        decision = "futility_stop",
        z1 = z1_obs,
        w1_ref = w1_ref,
        stage2_design = NULL,
        conditional_power = 0,
        recalibrated = recalibrate,
        message = "No designs feasible within budget"
      ), class = "interim_result"))
    }

    criterion <- ifelse(feasible, cp_vec, -Inf)
    best_idx <- which.max(criterion)
    best_cp <- cp_vec[best_idx]

    if (best_cp < 1e-6) {
      if (verbose) cat("\n** RECOMMENDATION: Stop for futility (negligible CP) **\n")
      return(structure(list(
        decision = "futility_stop",
        z1 = z1_obs,
        w1_ref = w1_ref,
        stage2_design = as.list(design_grid[best_idx, ]),
        conditional_power = best_cp,
        recalibrated = recalibrate,
        message = "Conditional power negligible within budget"
      ), class = "interim_result"))
    }

  } else {
    criterion <- cp_vec - lambda_use * cost_vec
    best_idx <- which.max(criterion)
    best_cp <- cp_vec[best_idx]
    best_criterion <- criterion[best_idx]

    if (verbose) cat(sprintf("Lambda = %.4e\n", lambda_use))

    if (best_criterion < 0) {
      if (verbose) cat("\n** RECOMMENDATION: Stop for futility **\n")
      return(structure(list(
        decision = "futility_stop",
        z1 = z1_obs,
        w1_ref = w1_ref,
        lambda = lambda_use,
        stage2_design = NULL,
        conditional_power = best_cp,
        criterion = best_criterion,
        recalibrated = recalibrate,
        message = sprintf("Best criterion = %.4f < 0", best_criterion)
      ), class = "interim_result"))
    }
  }

  # --- Step 6: Build recommendation ---
  best_design <- as.list(design_grid[best_idx, ])
  best_mod <- model_list[[best_idx]]

  if (verbose) {
    cat(sprintf("\n** RECOMMENDATION: Continue to stage 2 **\n"))
    cat(sprintf("Conditional power: %.1f%%\n", best_cp * 100))
    cat(sprintf("Stage 2 cost: %.1f\n", cost_vec[best_idx]))
    cat("Optimal stage 2 design:\n")
    for (nm in names(best_design)) {
      cat(sprintf("  %s = %s\n", nm, format(best_design[[nm]])))
    }
    cat(sprintf("I2_eff (updated): %.4f\n", I2_eff_vec[best_idx]))
    cat(sprintf("I2_eff (planned):  %.4f\n",
                r$results$models$list[[best_idx]]$I2_eff))
    if (recalibrate && length(theta_hat) > 0) {
      cat(sprintf("Lambda (planned):      %.4e\n", r$lambda %||% NA))
      cat(sprintf("Lambda (recalibrated): %.4e\n", lambda_use %||% NA))
    }
  }

  all_designs <- design_grid
  all_designs$cp <- cp_vec
  all_designs$cost <- cost_vec
  all_designs$I2_eff <- I2_eff_vec
  all_designs$selected <- seq_len(n_designs) == best_idx

  structure(list(
    decision = "continue",
    z1 = z1_obs,
    w1_ref = w1_ref,
    w2 = w2_ref,
    method = method,
    lambda = lambda_use,
    lambda_planned = if (method == "lambda") r$lambda else NA,
    cost_cap = cost_cap_use,
    stage2_design = best_design,
    conditional_power = best_cp,
    stage2_cost = cost_vec[best_idx],
    I2_eff_updated = I2_eff_vec[best_idx],
    theta_hat = theta_hat,
    delta = delta,
    recalibrated = recalibrate,
    all_designs = all_designs,
    planned = list(
      w1_ref = w1_ref,
      efficacy_boundary = efficacy_boundary,
      b1 = b1_planned
    ),
    message = "Continue to stage 2"
  ), class = "interim_result")
}


#' @export
print.interim_result <- function(x, ...) {
  cat("Interim Analysis Result\n")
  cat(sprintf("  Decision: %s\n", toupper(gsub("_", " ", x$decision))))
  cat(sprintf("  z1 observed: %.3f\n", x$z1))

  if (x$decision == "efficacy_stop") {
    cat(sprintf("  Weighted statistic: %.3f (critical value: %.3f)\n",
                x$weighted_statistic, x$critical_value))
  }

  if (x$decision == "continue") {
    cat(sprintf("  Conditional power: %.1f%%\n", x$conditional_power * 100))
    cat(sprintf("  Stage 2 cost: %.1f\n", x$stage2_cost))
    if (x$recalibrated && !is.na(x$lambda_planned)) {
      cat(sprintf("  Lambda: %.4e (planned: %.4e)\n", x$lambda, x$lambda_planned))
    }
    cat("  Stage 2 design:\n")
    for (nm in names(x$stage2_design)) {
      cat(sprintf("    %s = %s\n", nm, format(x$stage2_design[[nm]])))
    }
    if (length(x$theta_hat) > 0) {
      cat("  Updated parameters used:\n")
      for (nm in names(x$theta_hat)) {
        cat(sprintf("    %s = %s\n", nm, format(x$theta_hat[[nm]], digits = 4)))
      }
    }
  }

  if (x$decision == "futility_stop") {
    cat(sprintf("  Conditional power: %.1f%%\n", x$conditional_power * 100))
  }

  invisible(x)
}

#' Visualise interim decision rules under different parameter scenarios
#'
#' For each scenario of auxiliary parameters, computes the optimal stage 2
#' design across a grid of z1 values using pre-planned weights. Useful for
#' showing sensitivity of decision rules to interim estimates.
#'
#' @param planned An adaptive_results object (single design)
#' @param scenarios Named list of scenarios, each a named list of parameter values
#'   e.g. list("ICC = 0.03" = list(icc = 0.03), "ICC = 0.05" = list(icc = 0.05))
#' @param delta Treatment effect for CP (NULL = use planned)
#' @param z1_range Range of z1 values to evaluate
#' @param n_z1 Number of z1 points
#' @param type What to plot: "cost", "cp", "design", "all"
#' @param show_density Overlay density of z1 under H1
#'
#' @return A ggplot or list of ggplots
interim_sensitivity <- function(planned,
                                scenarios,
                                delta = NULL,
                                z1_range = c(-3, 4),
                                n_z1 = 200,
                                recalibrate = TRUE,
                                target_power = 0.8,
                                ...) {

  if (!inherits(planned, "adaptive_results")) {
    stop("planned must be an adaptive_results object")
  }
  if (is.null(scenarios)) stop("Provide at least one scenario in 'scenarios'")

  design <- planned$design
  spec <- design$spec
  rho <- design$rho

  r <- planned$raw
  w1_ref <- r$results$params$w1_ref
  w2_ref <- sqrt(1 - w1_ref^2)
  efficacy_boundary <- r$results$params$efficacy_boundary
  t_crit <- r$results$params$t_crit %||% 1.96
  b1_planned <- r$results$params$b1
  mu1_planned <- r$results$params$mu1
  method <- r$method %||% "lambda"

  if (is.null(delta)) delta <- b1_planned

  z1_grid <- seq(z1_range[1], z1_range[2], length.out = n_z1)
  dz <- diff(z1_grid[1:2])

  # Stage 1 params
  s1 <- design$stage1_grid[1, , drop = FALSE]
  design_grid <- design$stage2_grid_fn(s1)
  cost_fn <- generate_cost_fn(spec, "stage2")
  cost_fn_total <- generate_cost_fn(spec, "total")
  resource_vars <- spec$stage2_params

  # === Evaluate each scenario ===
  scenario_results <- lapply(names(scenarios), function(sc_name) {

    theta_hat <- scenarios[[sc_name]]

    # Build fixed params with updated estimates
    fixed_params <- design$fixed_params
    for (nm in names(s1)) fixed_params[[nm]] <- s1[[nm]]
    for (nm in names(theta_hat)) fixed_params[[nm]] <- theta_hat[[nm]]

    # Recompute models under this scenario
    model_fn <- generate_model_factory(spec)()

    if (recalibrate) {
      # Re-run full optimisation with locked weights
      recal <- find_lambda_for_power(
        target_power = target_power,
        design_grid = design_grid,
        model_fn = model_fn,
        fixed_params = fixed_params,
        cost_fn = cost_fn,
        cost_params = list(rho = rho),
        resource_vars = resource_vars,
        method = method,
        z1_range = z1_range,
        n_quad = n_z1,
        verbose = TRUE,
        w1_override = w1_ref,
        efficacy_override = efficacy_boundary,
        t_crit_override = t_crit,
        ...
      )

      # The results are already in plot_decision_rules-compatible format
      raw_compatible <- recal$results

      # Build decisions data frame from optimal_designs
      opt <- recal$results$quadrature$optimal_designs

      list(
        scenario = sc_name,
        theta = theta_hat,
        optimal_designs = opt,
        raw_compatible = raw_compatible,
        lambda = recal$lambda,
        cost_cap = recal$cost_cap,
        power = recal$power,
        converged = recal$converged
      )

    } else {
      # Original behaviour: use planned lambda, just recompute models
      model_list <- precompute_models(design_grid, model_fn, fixed_params,
                                      parallel = FALSE, n_cores = NULL)

      n_designs <- length(model_list)
      I2_vec <- sapply(model_list, function(d) d$I2_eff)
      cost_vec <- sapply(model_list, function(d) do.call(cost_fn, c(list(d), list(rho = rho))))

      I1_eff_scenario <- model_list[[1]]$I1_eff
      mu1_scenario <- delta * sqrt(I1_eff_scenario)

      mod1 <- model_list[[1]]
      cost_s1 <- cost_fn_total(mod1, rho) - cost_fn(mod1, rho)

      # CP matrix: n_z1 x n_designs
      cp_mat <- sapply(seq_len(n_designs), function(j) {
        ncp <- delta * sqrt(I2_vec[j])
        upper <- (t_crit - w1_ref * z1_grid) / w2_ref
        lower <- (-t_crit - w1_ref * z1_grid) / w2_ref
        pnorm(ncp - upper) + pnorm(lower - ncp)
      })

      # Apply decision criterion
      if (method == "cost_cap") {
        cost_cap <- r$cost_cap
        feasible <- cost_vec <= cost_cap
        cp_mat_feas <- cp_mat
        cp_mat_feas[, !feasible] <- -Inf
        best_idx <- max.col(cp_mat_feas)
        best_cp <- cp_mat[cbind(seq_along(z1_grid), best_idx)]
        futility <- best_cp < 1e-6
      } else {
        lambda <- r$lambda
        criterion_mat <- cp_mat - matrix(lambda * cost_vec,
                                         nrow = n_z1, ncol = n_designs,
                                         byrow = TRUE)
        best_idx <- max.col(criterion_mat)
        best_criterion <- criterion_mat[cbind(seq_along(z1_grid), best_idx)]
        best_cp <- cp_mat[cbind(seq_along(z1_grid), best_idx)]
        futility <- best_criterion < 0
      }

      efficacy <- abs(z1_grid) > efficacy_boundary
      continue <- !efficacy & !futility

      best_cost_s2 <- cost_vec[best_idx]
      best_cost_total <- cost_s1 + best_cost_s2

      opt <- data.frame(
        z1 = z1_grid,
        efficacy_stop = efficacy,
        futility_stop = futility,
        continue = continue,
        cp = ifelse(continue, best_cp, ifelse(efficacy, 1, 0)),
        cost = ifelse(continue, best_cost_total, cost_s1),
        design_idx = ifelse(continue, best_idx, NA)
      )

      for (col in names(design_grid)) {
        opt[[col]] <- ifelse(continue, design_grid[[col]][best_idx], NA)
      }
      for (var in resource_vars) {
        vals <- sapply(model_list, function(d) d[[var]] %||% NA)
        opt[[var]] <- ifelse(continue, vals[best_idx], NA)
      }

      # Quadrature weights and power
      weights <- dnorm(z1_grid, mean = mu1_scenario) * dz
      p_eff <- sum(weights[efficacy])
      p_fut <- sum(weights[futility])
      p_cont <- sum(weights[continue])
      power_s2 <- sum(weights[continue] * opt$cp[continue])
      upper_tail <- pnorm(max(z1_grid), mean = mu1_scenario, lower.tail = FALSE)
      total_power <- p_eff + power_s2 + upper_tail

      raw_compatible <- list(
        quadrature = list(
          optimal_designs = opt,
          weights = weights
        ),
        models = list(
          design_grid = design_grid
        ),
        params = list(
          mu1 = mu1_scenario,
          w1_ref = w1_ref,
          efficacy_boundary = efficacy_boundary,
          t_crit = t_crit,
          b1 = delta
        ),
        power = total_power,
        probabilities = list(
          efficacy_stop = p_eff,
          futility_stop = p_fut,
          continue = p_cont
        )
      )

      list(
        scenario = sc_name,
        theta = theta_hat,
        optimal_designs = opt,
        raw_compatible = raw_compatible,
        lambda = if (method == "lambda") r$lambda else NA,
        cost_cap = if (method == "cost_cap") r$cost_cap else NA,
        power = total_power,
        converged = TRUE
      )
    }
  })
  names(scenario_results) <- names(scenarios)

  # --- Build combined decisions data frame ---
  df <- do.call(rbind, lapply(scenario_results, function(sr) {
    opt <- sr$optimal_designs
    opt$scenario <- sr$scenario
    opt$outcome <- factor(
      ifelse(opt$efficacy_stop, "Efficacy stop",
             ifelse(opt$futility_stop, "Futility stop", "Continue")),
      levels = c("Efficacy stop", "Continue", "Futility stop")
    )
    opt
  }))
  df$scenario <- factor(df$scenario, levels = names(scenarios))

  structure(
    list(
      decisions = df,
      scenarios = scenarios,
      scenario_results = scenario_results,
      planned = list(
        w1_ref = w1_ref,
        w2_ref = w2_ref,
        mu1 = mu1_planned,
        efficacy_boundary = efficacy_boundary,
        t_crit = t_crit,
        delta = delta,
        method = method,
        lambda = if (method == "lambda") r$lambda else NA,
        cost_cap = if (method == "cost_cap") r$cost_cap else NA
      ),
      design_grid = design_grid,
      resource_vars = resource_vars,
      z1_grid = z1_grid,
      recalibrate = recalibrate,
      target_power = target_power
    ),
    class = "interim_sensitivity"
  )
}


#' @export
plot.interim_sensitivity <- function(x,
                                     type = c("decision", "cost", "cp",
                                              "design", "all"),
                                     show_density = TRUE,
                                     ...) {

  type <- match.arg(type)

  df <- x$decisions
  z1_grid <- x$z1_grid
  mu1 <- x$planned$mu1
  efficacy_boundary <- x$planned$efficacy_boundary

  # === Decision rules: delegate to plot_decision_rules per scenario ===
  if (type == "decision" || type == "all") {
    decision_plots <- lapply(names(x$scenarios), function(sc_name) {
      sr <- x$scenario_results[[sc_name]]
      fake <- sr$raw_compatible
      plot_decision_rules(fake, title = sc_name, show_density = show_density)
    })
    names(decision_plots) <- names(x$scenarios)

    if (requireNamespace("patchwork", quietly = TRUE) && length(decision_plots) > 1) {
      p_decision <- patchwork::wrap_plots(decision_plots, nrow = 1)
    } else {
      p_decision <- decision_plots
    }

    if (type == "decision") return(p_decision)
  }

  # === Overlay plots (scenarios on same axes) ===
  n_scenarios <- length(x$scenarios)
  palette_scenarios <- c("#264653", "#2a9d8f", "#e9c46a",
                         "#f4a261", "#e76f51")[1:min(n_scenarios, 5)]
  if (n_scenarios > 5) palette_scenarios <- scales::hue_pal()(n_scenarios)

  dens_df <- data.frame(z1 = z1_grid, density = dnorm(z1_grid - mu1))

  add_density <- function(p, y_max) {
    if (show_density) {
      dens_scaled <- dens_df$density / max(dens_df$density) * y_max * 0.25
      p + geom_ribbon(data = data.frame(z1 = z1_grid, ymin = 0,
                                        ymax = dens_scaled),
                      aes(x = z1, ymin = ymin, ymax = ymax),
                      fill = "grey85", alpha = 0.4, inherit.aes = FALSE)
    } else p
  }

  add_boundaries <- function(p) {
    p + geom_vline(xintercept = c(-efficacy_boundary, efficacy_boundary),
                   linetype = "dashed", colour = "grey50", linewidth = 0.5)
  }

  cont_df <- df[df$outcome == "Continue", ]

  # --- Cost ---
  p_cost <- {
    y_max <- max(cont_df$cost, na.rm = TRUE)
    p <- ggplot(df, aes(x = z1, y = cost, colour = scenario))
    p <- add_density(p, y_max)
    p <- p +
      geom_step(data = cont_df, linewidth = 0.8) +
      scale_colour_manual(values = palette_scenarios) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = expression(z[1]), y = "Total trial cost",
           colour = "Scenario") +
      theme_publication() +
      theme(legend.position = "bottom")
    add_boundaries(p)
  }

  # --- CP ---
  p_cp <- {
    p <- ggplot(df, aes(x = z1, y = cp, colour = scenario))
    p <- add_density(p, 1)
    p <- p +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = palette_scenarios) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, 1)) +
      labs(x = expression(z[1]), y = "Conditional power",
           colour = "Scenario") +
      theme_publication() +
      theme(legend.position = "bottom")
    add_boundaries(p)
  }

  # --- Design parameters ---
  design_cols <- intersect(names(x$design_grid), names(df))
  p_design_list <- lapply(design_cols, function(col) {
    if (all(is.na(cont_df[[col]]))) return(NULL)
    y_max <- max(cont_df[[col]], na.rm = TRUE)
    p <- ggplot(cont_df, aes(x = z1, y = .data[[col]], colour = scenario))
    p <- add_density(p, y_max)
    p <- p +
      geom_step(linewidth = 0.8) +
      scale_colour_manual(values = palette_scenarios) +
      labs(x = expression(z[1]), y = gsub("_", " ", col),
           colour = "Scenario") +
      theme_publication() +
      theme(legend.position = "bottom")
    add_boundaries(p)
  })
  p_design_list <- Filter(Negate(is.null), p_design_list)

  if (type == "cost") return(p_cost)
  if (type == "cp") return(p_cp)
  if (type == "design") {
    if (length(p_design_list) == 0) return(NULL)
    if (requireNamespace("patchwork", quietly = TRUE) && length(p_design_list) > 1) {
      return(
        patchwork::wrap_plots(p_design_list, ncol = 2) +
          patchwork::plot_layout(guides = "collect") &
          theme(legend.position = "bottom")
      )
    }
    return(p_design_list)
  }

  # type == "all"
  all_plots <- c(list(cost = p_cost, cp = p_cp), p_design_list)
  if (exists("p_decision")) all_plots <- c(list(decision = p_decision), all_plots)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(all_plots, ncol = 2) +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    all_plots
  }
}
