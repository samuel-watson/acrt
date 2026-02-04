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

  df[[group_var]] <- factor(df[[group_var]])
  n_groups <- length(unique(df[[group_var]]))

  palette_groups <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:min(n_groups, 5)]
  if (n_groups > 5) palette_groups <- scales::hue_pal()(n_groups)

  # Determine x variable
  x_var <- if ("power_stage1" %in% names(df)) {
    "power_stage1"
  } else if (length(stage1_vars) > 1) {
    stage1_vars[2]
  } else {
    group_var
  }

  if (!x_var %in% names(df)) x_var <- group_var

  # Cost ratio plot
  p_ratio <- ggplot(df, aes(x = .data[[x_var]])) +
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

  if (x_var == "power_stage1") {
    p_ratio <- p_ratio + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  }

  # Savings plot
  p_savings <- ggplot(df, aes(x = .data[[x_var]], y = E_savings_pct)) +
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

  if (x_var == "power_stage1") {
    p_savings <- p_savings + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  }

  # Frontier plot
  p_frontier <- ggplot(df, aes(x = E_cost, y = max_cost)) +
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
