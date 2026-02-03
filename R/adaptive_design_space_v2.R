# =============================================================================
# Design Space Exploration - Generalised Version
# =============================================================================

library(ggplot2)

# =============================================================================
# Theme for Publication-Quality Plots
# =============================================================================

theme_publication <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(size = rel(1), face = "plain", colour = "grey20"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10), angle = 90),
      axis.text = element_text(size = rel(0.85), colour = "grey30"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = rel(0.9), face = "bold", colour = "grey20"),
      legend.text = element_text(size = rel(0.8), colour = "grey30"),
      legend.position = "right",
      legend.key.size = unit(0.9, "lines"),
      legend.spacing.y = unit(0.3, "cm"),
      plot.title = element_text(size = rel(1.15), face = "bold", colour = "grey10",
                                margin = margin(b = 5), hjust = 0),
      plot.subtitle = element_text(size = rel(0.9), colour = "grey40",
                                   margin = margin(b = 15), hjust = 0),
      plot.caption = element_text(size = rel(0.75), colour = "grey50",
                                  hjust = 1, margin = margin(t = 10)),
      plot.margin = margin(15, 15, 10, 10)
    )
}

# =============================================================================
# Generic Design Space Exploration
# =============================================================================

#' Explore design space over stage 1 parameters
#'
#' @param stage1_grid Data frame of stage 1 parameter combinations to explore
#' @param target_power Target power for lambda search
#' @param design_grid_fn Function that returns stage 2 design grid
#' @param model_fn_factory Function that returns a model function (with fresh cache)
#' @param fixed_params_base Base fixed parameters (stage 1 params will be merged)
#' @param cost_fn Cost function
#' @param cost_params Fixed cost parameters (e.g., list(rho = 25))
#' @param resource_vars Resource variables to track
#' @param sample_size_fn Sample size function for computing E[N], N_max, etc.
#' @param summary_fn Optional function to compute additional summary metrics
#' @param verbose Print progress
#' @param ... Additional arguments passed to find_lambda_for_power
#'
#' @return List with summary data frame, individual results, and parameters
explore_design_space <- function(stage1_grid,
                                 target_power = 0.8,
                                 design_grid_fn,
                                 model_fn_factory,
                                 fixed_params_base,
                                 cost_fn,
                                 cost_params = list(rho = 25),
                                 resource_vars = NULL,
                                 sample_size_fn = NULL,
                                 summary_fn = NULL,
                                 verbose = TRUE,
                                 ...) {

  n_designs <- nrow(stage1_grid)

  if (verbose) cat(sprintf("Exploring %d stage 1 designs...\n\n", n_designs))

  # Storage for results
  results_list <- vector("list", n_designs)

  summary_list <- lapply(1:n_designs, function(i) {
    stage1_params <- stage1_grid[i, , drop = FALSE]

    if (verbose) {
      param_str <- paste(names(stage1_params), "=", stage1_params, collapse = ", ")
      cat(sprintf("=== Design %d/%d: %s ===\n", i, n_designs, param_str))
    }

    # Merge stage 1 params with fixed params
    fixed_params <- fixed_params_base
    for (nm in names(stage1_params)) {
      fixed_params[[nm]] <- stage1_params[[nm]]
    }

    # Get design grid for this stage 1 design
    design_grid <- design_grid_fn(stage1_params)

    # Create fresh model function (with fresh cache)
    model_fn <- model_fn_factory()

    # Find optimal lambda
    opt <- tryCatch({
      find_lambda_for_power(
        target_power = target_power,
        design_grid = design_grid,
        model_fn = model_fn,
        fixed_params = fixed_params,
        cost_fn = cost_fn,
        cost_params = cost_params,
        resource_vars = resource_vars,
        verbose = verbose,
        ...
      )
    }, error = function(e) {
      if (verbose) cat(sprintf("Error: %s\n", e$message))
      list(converged = FALSE, power = NA, lambda = NA,
           results = NULL, message = e$message)
    })

    results_list[[i]] <<- opt

    # Base summary: stage 1 params + convergence info
    base_summary <- as.data.frame(stage1_params)
    base_summary$converged <- opt$converged
    base_summary$lambda <- opt$lambda
    base_summary$power <- opt$power

    # Compute sample sizes if we have results and a sample_size_fn
    if (!is.null(opt$results) && !is.null(sample_size_fn)) {
      ss <- tryCatch({
        compute_sample_size(opt$results, sample_size_fn, hypothesis = "H1")
      }, error = function(e) {
        if (verbose) cat(sprintf("Sample size error: %s\n", e$message))
        NULL
      })

      if (!is.null(ss)) {
        # Add all summary metrics to base_summary
        for (j in seq_len(nrow(ss$summary))) {
          metric_name <- ss$summary$metric[j]
          # Convert metric names to valid R names with underscores
          col_name <- metric_name
          col_name <- gsub("\\[", "_", col_name)
          col_name <- gsub("\\]", "", col_name)
          col_name <- gsub("\\|", "_", col_name)
          col_name <- gsub("\\(", "_", col_name)
          col_name <- gsub("\\)", "", col_name)
          base_summary[[col_name]] <- ss$summary$value[j]
        }
      }
    }

    # Add expected resources from opt$results
    if (!is.null(opt$results) && !is.null(opt$results$expected)) {
      for (nm in names(opt$results$expected)) {
        base_summary[[paste0("E_", nm)]] <- opt$results$expected[[nm]]
      }
    }

    # Add probabilities
    if (!is.null(opt$results) && !is.null(opt$results$probabilities)) {
      base_summary$prob_efficacy <- opt$results$probabilities$efficacy_stop
      base_summary$prob_futility <- opt$results$probabilities$futility_stop
      base_summary$prob_continue <- opt$results$probabilities$continue
    }

    # Stage 1 power
    if (!is.null(opt$results) && !is.null(opt$results$params$mu1)) {
      mu1 <- opt$results$params$mu1
      base_summary$power_stage1 <- pnorm(-1.96 - mu1) + (1 - pnorm(1.96 - mu1))
    }

    # Call custom summary function if provided
    if (!is.null(summary_fn) && !is.null(opt$results)) {
      custom_metrics <- tryCatch({
        summary_fn(opt, stage1_params, design_grid)
      }, error = function(e) {
        if (verbose) cat(sprintf("Summary function error: %s\n", e$message))
        NULL
      })

      if (!is.null(custom_metrics)) {
        for (nm in names(custom_metrics)) {
          base_summary[[nm]] <- custom_metrics[[nm]]
        }
      }
    }

    if (verbose) cat("\n")

    base_summary
  })

  # FIX: Ensure all data frames have the same columns before rbind
  all_cols <- unique(unlist(lapply(summary_list, names)))

  summary_list <- lapply(summary_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    # Reorder to match all_cols
    df[, all_cols, drop = FALSE]
  })

  summary_df <- do.call(rbind, summary_list)

  list(
    summary = summary_df,
    results = results_list,
    params = list(
      stage1_grid = stage1_grid,
      target_power = target_power,
      cost_params = cost_params
    )
  )
}

# =============================================================================
# Convenience Wrappers for Common Cases
# =============================================================================

#' Create a standard (k1, m1) stage 1 grid
#'
#' @param k1_values Vector of k1 values (clusters per arm)
#' @param m1_values Vector of m1 values (individuals per cluster-period)
#' @return Data frame suitable for explore_design_space
make_stage1_grid <- function(k1_values, m1_values) {
  expand.grid(k1 = k1_values, m1 = m1_values)
}

#' Create a standard (m2, k2) stage 2 design grid function
#'
#' @param m2_range Range or sequence for m2
#' @param k2_range Range or sequence for k2
#' @param m2_fn Optional function to compute m2 range from stage 1 params
#' @return Function suitable for design_grid_fn argument
make_stage2_grid_fn <- function(m2_range = NULL, k2_range = 0:3, m2_fn = NULL) {
  function(stage1_params) {
    if (!is.null(m2_fn)) {
      m2_seq <- m2_fn(stage1_params)
    } else if (!is.null(m2_range)) {
      m2_seq <- m2_range
    } else {
      # Default: 0.5*m1 to 2.5*m1
      m1 <- stage1_params$m1
      m2_seq <- seq(max(10, floor(0.5 * m1)), ceiling(2.5 * m1), by = 5)
    }
    expand.grid(m2 = m2_seq, k2 = k2_range)
  }
}

# =============================================================================
# Pareto Frontier Functions
# =============================================================================

#' Find constrained Pareto frontier
#'
#' @param df Data frame (typically exploration$summary or comparison$comparison)
#' @param objectives Named list of objectives with directions.
#'        E.g., list(E_N = "min", N_max = "min") or c(E_N = "min", prob_efficacy = "max")
#' @param constraints Named list of constraints. Each element is either:
#'        - A numeric vector of length 2 for range constraints: c(min, max)
#'        - A vector of allowed values for set membership
#' @param require_converged Only include converged designs (default TRUE)
#' @param preset Optional preset: "cost_efficient", "cluster_efficient", "sample_efficient"
#' @return Data frame of Pareto-optimal designs with attributes
find_constrained_pareto <- function(df,
                                    objectives = list(E_N = "min", N_max = "min"),
                                    constraints = NULL,
                                    require_converged = TRUE,
                                    preset = NULL) {

  # Handle presets
  if (!is.null(preset)) {
    objectives <- switch(preset,
      cost_efficient = list(E_cost = "min", max_cost = "min"),
      cluster_efficient = list(E_N = "min", E_total_clusters = "min"),
      sample_efficient = list(E_N = "min", N_max = "min"),
      stop("Unknown preset: ", preset)
    )
  }

  # Convert named character vector to list if needed
  if (!is.list(objectives)) {
    objectives <- as.list(objectives)
  }

  obj_names <- names(objectives)
  obj_directions <- unlist(objectives)

  # Validate objectives exist in df
  missing <- setdiff(obj_names, names(df))
  if (length(missing) > 0) {
    stop("Objectives not found in data: ", paste(missing, collapse = ", "))
  }

  # Filter for converged designs
  if (require_converged && "converged" %in% names(df)) {
    df <- df[df$converged, ]
  }

  # Apply constraints
  if (!is.null(constraints)) {
    for (var in names(constraints)) {
      if (!var %in% names(df)) {
        warning("Constraint variable not found: ", var)
        next
      }

      constraint <- constraints[[var]]

      if (length(constraint) == 2 && is.numeric(constraint)) {
        # Range constraint
        df <- df[df[[var]] >= constraint[1] & df[[var]] <= constraint[2], ]
      } else {
        # Set membership
        df <- df[df[[var]] %in% constraint, ]
      }
    }
  }

  # Remove NA in objective columns
  complete_idx <- complete.cases(df[, obj_names, drop = FALSE])
  df <- df[complete_idx, ]

  if (nrow(df) == 0) {
    warning("No designs remain after filtering")
    result <- df
    attr(result, "objectives") <- objectives
    attr(result, "n_dominated") <- 0
    return(result)
  }

  if (nrow(df) == 1) {
    result <- df
    attr(result, "objectives") <- objectives
    attr(result, "n_dominated") <- 0
    return(result)
  }

  # Convert to minimisation problem
  obj_matrix <- as.matrix(df[, obj_names, drop = FALSE])
  for (i in seq_along(obj_names)) {
    if (obj_directions[i] == "max") {
      obj_matrix[, i] <- -obj_matrix[, i]
    }
  }

  # Find Pareto-optimal designs
  is_dominated <- function(i) {
    obj_i <- obj_matrix[i, ]
    for (j in seq_len(nrow(df))) {
      if (i == j) next
      obj_j <- obj_matrix[j, ]
      # j dominates i if j <= i on all and j < i on at least one
      if (all(obj_j <= obj_i) && any(obj_j < obj_i)) return(TRUE)
    }
    FALSE
  }

  dominated <- sapply(seq_len(nrow(df)), is_dominated)
  n_dominated <- sum(dominated)

  result <- df[!dominated, ]
  result <- result[order(result[[obj_names[1]]]), ]

  attr(result, "objectives") <- objectives
  attr(result, "n_dominated") <- n_dominated

  result
}

#' Compute Pareto ranks for all designs
#'
#' @param df Data frame of designs
#' @param objectives Named list of objectives with directions
#' @param max_rank Maximum rank to compute (NULL for all)
#' @return Data frame with pareto_rank column added
compute_pareto_ranks <- function(df, objectives, max_rank = NULL) {

  obj_names <- names(objectives)
  obj_directions <- unlist(objectives)

  # Convert to minimisation
  obj_matrix <- as.matrix(df[, obj_names, drop = FALSE])
  for (i in seq_along(obj_names)) {
    if (obj_directions[i] == "max") {
      obj_matrix[, i] <- -obj_matrix[, i]
    }
  }

  # Iteratively find Pareto frontiers
  remaining_idx <- seq_len(nrow(df))
  df$pareto_rank <- NA_integer_
  rank <- 0

  while (length(remaining_idx) > 0) {
    rank <- rank + 1

    if (!is.null(max_rank) && rank > max_rank) break

    # Find non-dominated among remaining
    is_dominated_among_remaining <- function(i) {
      obj_i <- obj_matrix[i, ]
      for (j in remaining_idx) {
        if (i == j) next
        obj_j <- obj_matrix[j, ]
        if (all(obj_j <= obj_i) && any(obj_j < obj_i)) return(TRUE)
      }
      FALSE
    }

    dominated <- sapply(remaining_idx, is_dominated_among_remaining)
    frontier_idx <- remaining_idx[!dominated]

    df$pareto_rank[frontier_idx] <- rank
    remaining_idx <- remaining_idx[dominated]
  }

  df[order(df$pareto_rank, df[[obj_names[1]]]), ]
}

# =============================================================================
# Plotting Functions
# =============================================================================

#' Plot Pareto frontier
#'
#' @param pareto_df Output from find_constrained_pareto
#' @param all_designs Optional data frame of all designs (for showing dominated points)
#' @param label_var Variable to use for point labels
#' @param show_frontier_line Show step line connecting frontier points
#' @return ggplot object
plot_pareto <- function(pareto_df,
                        all_designs = NULL,
                        label_var = NULL,
                        show_frontier_line = TRUE) {

  obj_list <- attr(pareto_df, "objectives")

  if (is.null(obj_list)) {
    stop("No 'objectives' attribute found. Ensure pareto_df comes from find_constrained_pareto().")
  }

  if (length(obj_list) < 2) {
    stop(sprintf("Need at least 2 objectives for plotting, found %d", length(obj_list)))
  }

  obj_names <- names(obj_list)[1:2]
  obj_directions <- unlist(obj_list)[1:2]

  x_var <- obj_names[1]
  y_var <- obj_names[2]

  p <- ggplot()

  # Show dominated designs if provided
  if (!is.null(all_designs)) {
    if (nrow(pareto_df) > 0 && all(c(x_var, y_var) %in% names(all_designs))) {
      all_designs$.on_frontier <- interaction(all_designs[[x_var]], all_designs[[y_var]]) %in%
        interaction(pareto_df[[x_var]], pareto_df[[y_var]])

      dominated_designs <- all_designs[!all_designs$.on_frontier, ]

      if (nrow(dominated_designs) > 0) {
        p <- p +
          geom_point(data = dominated_designs,
                     aes(x = .data[[x_var]], y = .data[[y_var]]),
                     colour = "grey70", size = 2, alpha = 0.5)
      }
    }
  }

  # Pareto frontier line
  if (show_frontier_line && nrow(pareto_df) > 1) {
    pareto_sorted <- pareto_df[order(pareto_df[[x_var]]), ]

    p <- p +
      geom_step(data = pareto_sorted,
                aes(x = .data[[x_var]], y = .data[[y_var]]),
                colour = "#264653", linewidth = 0.8, alpha = 0.6,
                direction = "vh")
  }

  # Pareto points
  p <- p +
    geom_point(data = pareto_df,
               aes(x = .data[[x_var]], y = .data[[y_var]]),
               fill = "#2a9d8f", colour = "white",
               shape = 21, size = 4, stroke = 0.8)

  # Labels if requested
  if (!is.null(label_var) && label_var %in% names(pareto_df)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = pareto_df,
          aes(x = .data[[x_var]], y = .data[[y_var]],
              label = .data[[label_var]]),
          size = 3, colour = "grey30",
          segment.colour = "grey70", segment.size = 0.3
        )
    } else {
      p <- p +
        geom_text(data = pareto_df,
                  aes(x = .data[[x_var]], y = .data[[y_var]],
                      label = .data[[label_var]]),
                  hjust = -0.2, vjust = -0.2, size = 3, colour = "grey30")
    }
  }

  # Axis labels with direction indicators
  make_axis_label <- function(var_name, direction) {
    display_name <- gsub("_", " ", var_name)
    dir_label <- if (direction == "min") "minimise" else "maximise"
    sprintf("%s (%s)", display_name, dir_label)
  }

  x_label <- make_axis_label(x_var, obj_directions[1])
  y_label <- make_axis_label(y_var, obj_directions[2])

  n_designs <- nrow(pareto_df)
  subtitle <- if (n_designs == 1) {
    "Single optimal design"
  } else {
    sprintf("%d non-dominated designs", n_designs)
  }

  p <- p +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x = x_label,
      y = y_label,
      title = "Pareto frontier",
      subtitle = subtitle
    ) +
    theme_publication()

  p
}

#' Plot design space exploration results
#'
#' @param exploration Output from explore_design_space
#' @param x_var Variable for x-axis (default: first column of stage1_grid)
#' @param y_var Variable for y-axis (default: "E_N")
#' @param colour_var Variable for colour grouping (default: second column of stage1_grid)
#' @param label_var Variable for point labels (optional)
#' @param facet_var Variable for faceting (optional)
plot_exploration <- function(exploration,
                             x_var = NULL,
                             y_var = "E_N",
                             colour_var = NULL,
                             label_var = NULL,
                             facet_var = NULL) {

  df <- exploration$summary
  stage1_vars <- names(exploration$params$stage1_grid)

  # Default x and colour from stage1_grid
  if (is.null(x_var)) x_var <- stage1_vars[1]
  if (is.null(colour_var) && length(stage1_vars) > 1) colour_var <- stage1_vars[2]

  # Ensure colour_var is a factor
  if (!is.null(colour_var)) {
    df[[colour_var]] <- factor(df[[colour_var]])
    n_colours <- length(unique(df[[colour_var]]))
    palette <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:min(n_colours, 5)]
    if (n_colours > 5) {
      palette <- viridis::viridis(n_colours)
    }
  }

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]]))

  if (!is.null(colour_var)) {
    p <- p +
      geom_line(aes(colour = .data[[colour_var]], group = .data[[colour_var]]),
                linewidth = 0.9) +
      geom_point(aes(colour = .data[[colour_var]], shape = .data[[colour_var]]),
                 size = 3, fill = "white", stroke = 0.8) +
      scale_colour_manual(values = palette) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)[1:n_colours])
  } else {
    p <- p +
      geom_line(linewidth = 0.9, colour = "#264653") +
      geom_point(size = 3, colour = "#264653")
  }

  if (!is.null(label_var) && label_var %in% names(df)) {
    p <- p +
      geom_text(aes(label = .data[[label_var]]),
                hjust = -0.4, vjust = -0.4, size = 2.8, colour = "grey40")
  }

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }

  # Format y-axis based on variable type
  y_label <- gsub("_", " ", y_var)
  if (grepl("^prob|^P_", y_var)) {
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  } else if (grepl("^E_|^N_", y_var) || y_var %in% c("E_N", "N_max")) {
    p <- p + scale_y_continuous(labels = scales::comma)
  }

  p <- p +
    labs(
      x = gsub("_", " ", x_var),
      y = y_label,
      title = sprintf("%s by design parameters", y_label),
      subtitle = sprintf("Target power = %g%%", exploration$params$target_power * 100)
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  p
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
