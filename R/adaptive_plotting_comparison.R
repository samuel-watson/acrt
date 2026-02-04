# =============================================================================
# Plotting and Comparison Functions
# =============================================================================

# =============================================================================
# Decision Rules Plotting
# =============================================================================

#' Plot decision rules across z1 values
#'
#' Creates publication-quality visualisation of adaptive trial decision rules.
#'
#' @param results Output from adaptive_analysis or find_lambda_for_power
#' @param design_vars Which design variables to show (default: all)
#' @param design_index For explored results, which design to plot (default: 1)
#' @param show_density Show z1 density under H1 (default TRUE)
#' @param title Custom title
#' @return A ggplot or patchwork object
plot_decision_rules <- function(results,
                                design_vars = NULL,
                                design_index = 1,
                                show_density = TRUE,
                                title = NULL) {



  # === Extract results from various input formats ===
  raw <- extract_raw_results(results, design_index)

  if (is.null(raw) || is.null(raw$quadrature) || is.null(raw$quadrature$optimal_designs)) {
    stop("No optimal designs found in results.")
  }

  opt <- raw$quadrature$optimal_designs
  design_grid <- raw$models$design_grid
  mu1 <- raw$params$mu1 %||% 0
  weights <- raw$quadrature$weights

  # Get design variables
  if (is.null(design_vars)) {
    design_vars <- setdiff(names(design_grid), c("row_id"))
  }

  # === Build plot data ===
  plot_data <- opt
  plot_data$outcome <- factor(
    ifelse(plot_data$efficacy_stop, "Efficacy",
           ifelse(plot_data$futility_stop, "Futility", "Continue")),
    levels = c("Efficacy", "Futility", "Continue")
  )

  # Calculate density under H1
  plot_data$density_h1 <- dnorm(plot_data$z1, mean = mu1, sd = 1)

  # P(reject | z1)
  plot_data$p_reject <- ifelse(
    plot_data$efficacy_stop, 1,
    ifelse(plot_data$futility_stop, 0, plot_data$cp)
  )

  # === Find decision boundaries ===
  boundaries <- find_decision_boundaries(plot_data)

  # === Calculate summary probabilities ===
  if (!is.null(weights) && length(weights) == nrow(plot_data)) {
    p_efficacy <- sum(weights[plot_data$outcome == "Efficacy"])
    p_futility <- sum(weights[plot_data$outcome == "Futility"])
    p_continue <- sum(weights[plot_data$outcome == "Continue"])
  } else {
    # Approximate from density
    total_dens <- sum(plot_data$density_h1)
    p_efficacy <- sum(plot_data$density_h1[plot_data$outcome == "Efficacy"]) / total_dens
    p_futility <- sum(plot_data$density_h1[plot_data$outcome == "Futility"]) / total_dens
    p_continue <- sum(plot_data$density_h1[plot_data$outcome == "Continue"]) / total_dens
  }

  # Overall power
  # === Get stored power and probabilities ===
  total_power <- raw$power
  if (is.list(total_power)) {
    total_power <- total_power$total %||% total_power[[1]]
  }

  probs <- raw$probabilities
  p_efficacy <- probs$efficacy_stop
  if (length(p_efficacy) > 1) p_efficacy <- p_efficacy[1]
  p_futility <- probs$futility_stop
  if (length(p_futility) > 1) p_futility <- p_futility[1]
  p_continue <- probs$continue
  if (length(p_continue) > 1) p_continue <- p_continue[1]

  # Get efficacy boundary info for display
  w1_ref <- raw$params$w1_ref %||% NULL
  efficacy_boundary <- raw$params$efficacy_boundary %||% NULL

  # === Colour palette ===
  palette <- c(
    "Efficacy"  = "#1a759f",
    "Futility"  = "#d62828",
    "Continue"  = "#2a9d8f"
  )

  # === Main plot ===
  p_main <- plot_decisions_main(plot_data, boundaries, palette, mu1,
                                p_efficacy, p_futility, p_continue,
                                total_power, show_density, title)

  # === Design variable plots ===
  design_plots <- create_design_plots(plot_data, design_vars, boundaries, palette)

  # === Combine ===
  if (length(design_plots) == 0) {
    return(p_main)
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p_main
    for (dp in design_plots) {
      combined <- combined / dp
    }

    combined +
      patchwork::plot_layout(heights = c(2.5, rep(1, length(design_plots)))) +
      patchwork::plot_annotation(
        theme = theme(plot.margin = margin(5, 5, 5, 5))
      )
  } else {
    list(main = p_main, design = design_plots)
  }
}

#' Extract raw results from various input formats
#' @keywords internal
extract_raw_results <- function(results, design_index = 1) {
  if (inherits(results, "adaptive_results")) {
    if (results$explored) {
      results$raw$results[[design_index]]$results
    } else {
      results$raw$results
    }
  } else if (!is.null(results$quadrature)) {
    results
  } else if (!is.null(results$results)) {
    results$results
  } else {
    NULL
  }
}

#' Find decision boundaries from plot data
#' @keywords internal
find_decision_boundaries <- function(plot_data) {
  boundaries <- list()

  for (outcome in c("Efficacy", "Futility", "Continue")) {
    z1_vals <- plot_data$z1[plot_data$outcome == outcome]
    if (length(z1_vals) > 0) {
      prefix <- tolower(outcome)
      boundaries[[paste0(prefix, "_min")]] <- min(z1_vals, na.rm = TRUE)
      boundaries[[paste0(prefix, "_max")]] <- max(z1_vals, na.rm = TRUE)
    }
  }

  boundaries
}

#' Main decision plot showing regions and density
#' @keywords internal
plot_decisions_main <- function(plot_data, boundaries, palette, mu1,
                                p_efficacy, p_futility, p_continue,
                                total_power, show_density, title,
                                w1_ref = NULL, efficacy_boundary = NULL) {

  # Scale density to 0-1 range
  plot_data$density_scaled <- plot_data$density_h1 / max(plot_data$density_h1)

  # Sort by z1 for finding contiguous blocks
  plot_data <- plot_data[order(plot_data$z1), ]

  # Find contiguous blocks of each outcome
  blocks <- find_contiguous_blocks(plot_data$z1, plot_data$outcome)

  p <- ggplot(plot_data, aes(x = z1))

  # Background regions for each contiguous block
  for (i in seq_len(nrow(blocks))) {
    p <- p + annotate("rect",
                      xmin = blocks$z1_min[i],
                      xmax = blocks$z1_max[i],
                      ymin = -Inf, ymax = Inf,
                      fill = palette[blocks$outcome[i]], alpha = 0.12)
  }

  # Density curve
  if (show_density) {
    p <- p +
      geom_ribbon(aes(ymin = 0, ymax = density_scaled),
                  fill = "#264653", alpha = 0.2) +
      geom_line(aes(y = density_scaled), colour = "#264653", linewidth = 0.6)
  }

  # Decision outcome as colored bar at bottom
  p <- p +
    geom_tile(aes(y = -0.08, fill = outcome), height = 0.06, alpha = 0.9) +
    scale_fill_manual(values = palette, name = "Decision")

  # Region labels - one per block
  for (i in seq_len(nrow(blocks))) {
    mid_z1 <- (blocks$z1_min[i] + blocks$z1_max[i]) / 2
    outcome <- blocks$outcome[i]

    # Calculate probability for this specific block
    block_data <- plot_data[plot_data$z1 >= blocks$z1_min[i] &
                              plot_data$z1 <= blocks$z1_max[i] &
                              plot_data$outcome == outcome, ]

    # Only label if block is wide enough
    block_width <- blocks$z1_max[i] - blocks$z1_min[i]
    total_width <- max(plot_data$z1) - min(plot_data$z1)

    if (block_width / total_width > 0.08) {  # Only label if > 8% of total width
      prob <- switch(outcome,
                     "Efficacy" = p_efficacy,
                     "Futility" = p_futility,
                     "Continue" = p_continue)

      # If there are multiple blocks of same outcome, don't repeat percentage
      n_blocks_same <- sum(blocks$outcome == outcome)
      if (n_blocks_same > 1) {
        label_text <- outcome
      } else {
        label_text <- sprintf("%s\n%.0f%%", outcome, prob * 100)
      }

      p <- p + annotate("label",
                        x = mid_z1, y = 0.92,
                        label = label_text,
                        fill = palette[outcome], colour = "white",
                        fontface = "bold", size = 3,
                        label.padding = unit(0.3, "lines"))
    }
  }

  # Add overall probabilities to subtitle instead
  # Build subtitle
  subtitle <- sprintf(
    "Overall power: %.1f%%   |   μ₁ = %.2f",
    total_power * 100, mu1
  )

  if (!is.null(efficacy_boundary)) {
    subtitle <- paste0(subtitle, sprintf("   |   Efficacy: |z₁| > %.2f", efficacy_boundary))
  }

  p +
    scale_y_continuous(
      limits = c(-0.15, 1.05),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      x = expression(z[1]~"(stage 1 statistic)"),
      y = if (show_density) expression("Density of "*z[1]*" under "*H[1]) else NULL,
      title = title %||% "Interim Decision Rules",
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey30", size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.title.y = element_text(colour = "#264653"),
      plot.margin = margin(10, 15, 10, 10)
    )
}

#' Find contiguous blocks of outcomes with extended boundaries
#' @keywords internal
find_contiguous_blocks <- function(z1, outcome) {
  n <- length(z1)
  if (n == 0) return(data.frame(outcome = character(), z1_min = numeric(), z1_max = numeric()))

  # Identify where outcome changes
  outcome_char <- as.character(outcome)
  changes <- c(TRUE, outcome_char[-1] != outcome_char[-n])
  block_id <- cumsum(changes)

  n_blocks <- max(block_id)

  # Summarise each block
  blocks <- data.frame(
    block = 1:n_blocks
  )

  blocks$outcome <- sapply(blocks$block, function(b) outcome_char[block_id == b][1])
  blocks$z1_min_raw <- sapply(blocks$block, function(b) min(z1[block_id == b]))
  blocks$z1_max_raw <- sapply(blocks$block, function(b) max(z1[block_id == b]))

  # Extend boundaries to midpoints between blocks
  blocks$z1_min <- blocks$z1_min_raw
  blocks$z1_max <- blocks$z1_max_raw

  for (i in 1:n_blocks) {
    if (i > 1) {
      # Extend left boundary to midpoint with previous block
      blocks$z1_min[i] <- (blocks$z1_max_raw[i - 1] + blocks$z1_min_raw[i]) / 2
    }
    if (i < n_blocks) {
      # Extend right boundary to midpoint with next block
      blocks$z1_max[i] <- (blocks$z1_max_raw[i] + blocks$z1_min_raw[i + 1]) / 2
    }
  }

  blocks[, c("outcome", "z1_min", "z1_max")]
}

#' Create design variable plots
#' @keywords internal
create_design_plots <- function(plot_data, design_vars, boundaries, palette) {

  design_plots <- lapply(design_vars, function(var) {
    if (!var %in% names(plot_data)) return(NULL)

    df_var <- plot_data[plot_data$outcome == "Continue", ]
    if (nrow(df_var) == 0) return(NULL)
    if (all(is.na(df_var[[var]]))) return(NULL)

    var_label <- tools::toTitleCase(gsub("_", " ", var))

    p <- ggplot(df_var, aes(x = z1, y = .data[[var]])) +
      geom_ribbon(aes(ymin = min(.data[[var]], na.rm = TRUE),
                      ymax = .data[[var]]),
                  fill = "#2a9d8f", alpha = 0.15) +
      geom_step(linewidth = 0.9, colour = "#2a9d8f") +
      geom_point(size = 1.8, colour = "#2a9d8f", alpha = 0.7) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      scale_y_continuous(expand = expansion(mult = 0.1)) +
      labs(x = expression(z[1]), y = var_label) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 15, 5, 10)
      )

    p
  })

  Filter(Negate(is.null), design_plots)
}

#' Compact decision rules plot
#' @param results Output from adaptive_analysis or find_lambda_for_power
#' @param design_var Single design variable to show (default: first available)
#' @param design_index For explored results, which design to plot
#' @param title Custom title
#' @return A ggplot object
plot_decision_rules_compact <- function(results,
                                        design_var = NULL,
                                        design_index = 1,
                                        title = NULL) {



  raw <- extract_raw_results(results, design_index)

  if (is.null(raw) || is.null(raw$quadrature) || is.null(raw$quadrature$optimal_designs)) {
    stop("No optimal designs found in results.")
  }

  opt <- raw$quadrature$optimal_designs
  design_grid <- raw$models$design_grid
  mu1 <- raw$params$mu1 %||% 0

  # Get design variable
  if (is.null(design_var)) {
    available_vars <- setdiff(names(design_grid), c("row_id"))
    if (length(available_vars) == 0) {
      available_vars <- setdiff(names(opt), c("z1", "continue", "efficacy_stop",
                                              "futility_stop", "cp", "criterion"))
    }
    if (length(available_vars) > 0) design_var <- available_vars[1]
  }

  # Build data
  plot_data <- opt
  plot_data$outcome <- factor(
    ifelse(plot_data$efficacy_stop, "Efficacy",
           ifelse(plot_data$futility_stop, "Futility", "Continue")),
    levels = c("Efficacy", "Futility", "Continue")
  )

  plot_data$p_reject <- ifelse(
    plot_data$efficacy_stop, 1,
    ifelse(plot_data$futility_stop, 0, plot_data$cp)
  )

  # Design variable
  has_design_var <- !is.null(design_var) && design_var %in% names(plot_data)
  if (has_design_var) {
    plot_data$design_val <- plot_data[[design_var]]
  }

  palette <- c(
    "Efficacy"  = "#1a759f",
    "Futility"  = "#d62828",
    "Continue"  = "#2a9d8f"
  )

  boundaries <- find_decision_boundaries(plot_data)

  # Design variable range
  continue_data <- plot_data[plot_data$outcome == "Continue", ]
  if (has_design_var && nrow(continue_data) > 0 && !all(is.na(continue_data$design_val))) {
    design_range <- range(continue_data$design_val, na.rm = TRUE)
    if (design_range[1] == design_range[2]) design_range <- design_range + c(-0.5, 0.5)
  } else {
    design_range <- c(0, 1)
    has_design_var <- FALSE
  }

  z_range <- range(plot_data$z1)
  z_buffer <- diff(z_range) * 0.02

  p <- ggplot(plot_data, aes(x = z1))

  # Background regions
  if (!is.null(boundaries$futility_min)) {
    p <- p + annotate("rect",
                      xmin = boundaries$futility_min - z_buffer,
                      xmax = boundaries$futility_max + z_buffer,
                      ymin = -Inf, ymax = Inf,
                      fill = palette["Futility"], alpha = 0.1)
  }
  if (!is.null(boundaries$efficacy_min)) {
    p <- p + annotate("rect",
                      xmin = boundaries$efficacy_min - z_buffer,
                      xmax = boundaries$efficacy_max + z_buffer,
                      ymin = -Inf, ymax = Inf,
                      fill = palette["Efficacy"], alpha = 0.1)
  }

  # P(reject | z1)
  p <- p +
    geom_line(aes(y = p_reject), linewidth = 1, colour = "#264653") +
    geom_point(aes(y = p_reject, colour = outcome), size = 2.5) +
    scale_colour_manual(values = palette, name = "Decision")

  # Design variable
  if (has_design_var && nrow(continue_data) > 0) {
    continue_data$design_scaled <- scales::rescale(
      continue_data$design_val, to = c(0.1, 0.7), from = design_range
    )

    p <- p +
      geom_step(data = continue_data, aes(y = design_scaled),
                linewidth = 0.9, colour = "#e9c46a") +
      geom_point(data = continue_data, aes(y = design_scaled),
                 size = 1.5, colour = "#e9c46a", alpha = 0.7)
  }

  # Reference
  p <- p + geom_hline(yintercept = 0.8, linetype = "dotted", colour = "grey60")

  # Scales
  if (has_design_var) {
    var_label <- tools::toTitleCase(gsub("_", " ", design_var))
    p <- p + scale_y_continuous(
      name = expression(P(reject~H[0]~"|"~z[1])),
      labels = scales::percent,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(
        ~ scales::rescale(., from = c(0.1, 0.7), to = design_range),
        name = var_label,
        breaks = scales::pretty_breaks(n = 5)
      )
    )
  } else {
    p <- p + scale_y_continuous(
      name = expression(P(reject~H[0]~"|"~z[1])),
      labels = scales::percent,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    )
  }

  p +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      x = expression(z[1]~"(stage 1 statistic)"),
      title = title %||% "Adaptive Decision Rules"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      axis.title.y.left = element_text(colour = "#264653", size = 10),
      axis.title.y.right = element_text(colour = "#e9c46a", face = "bold", size = 10),
      axis.text.y.right = element_text(colour = "#e9c46a"),
      plot.title = element_text(face = "bold", size = 12),
      plot.margin = margin(10, 15, 10, 10)
    )
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# =============================================================================
# Sample Size Plotting
# =============================================================================

#' Plot sample size distribution
#' @param ss Output from compute_sample_size or compare_sample_sizes
#' @param type "histogram" (default), "cumulative", or "both"
#' @param show_outcomes Colour by outcome type (default TRUE)
plot_sample_size <- function(ss, type = c("histogram", "cumulative", "both"),
                             show_outcomes = TRUE) {

  type <- match.arg(type)


  # Handle compare_sample_sizes output
  if ("comparison" %in% names(ss)) {
    return(plot_sample_size_comparison(ss, type, show_outcomes))
  }

  # Colour palette
  palette_outcome <- c(
    "Efficacy stop" = "#2d6a4f",
    "Futility stop" = "#ae2012",
    "Continue"      = "#457b9d"
  )
  col_reference <- "grey60"

  dist <- ss$distribution

  dist$outcome <- factor(
    ifelse(dist$efficacy_stop, "Efficacy stop",
           ifelse(dist$futility_stop, "Futility stop", "Continue")),
    levels = c("Efficacy stop", "Futility stop", "Continue")
  )

  # Aggregate by unique sample sizes
  dist_agg <- aggregate(weight ~ n_total + outcome, data = dist, FUN = sum)

  # Summary statistics
  E_N <- ss$summary$value[ss$summary$metric == "E[N]"]
  N_max <- max(dist$n_total)

  # Histogram
  bar_width <- max(diff(range(dist_agg$n_total)) / 40, 1)

  p_hist <- ggplot(dist_agg, aes(x = n_total, y = weight)) +
    geom_vline(xintercept = E_N, linetype = "dashed",
               colour = col_reference, linewidth = 0.6) +
    annotate("text", x = E_N, y = Inf,
             label = sprintf("E[N] = %s", scales::comma(round(E_N))),
             hjust = -0.1, vjust = 2, size = 3, colour = "grey40", fontface = "italic")

  if (show_outcomes) {
    p_hist <- p_hist +
      geom_col(aes(fill = outcome), width = bar_width, alpha = 0.85) +
      scale_fill_manual(values = palette_outcome, name = "Outcome")
  } else {
    p_hist <- p_hist +
      geom_col(fill = "#457b9d", width = bar_width, alpha = 0.85)
  }

  p_hist <- p_hist +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      x = "Total sample size (N)",
      y = "Probability",
      title = "Sample size distribution",
      subtitle = sprintf("Under %s", ss$hypothesis)
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  # Cumulative plot
  dist_sorted <- dist[order(dist$n_total), ]
  dist_sorted$cum_prob <- cumsum(dist_sorted$weight) / sum(dist_sorted$weight)

  p_cdf <- ggplot(dist_sorted, aes(x = n_total, y = cum_prob)) +
    geom_step(linewidth = 0.9, colour = "#264653") +
    geom_hline(yintercept = c(0.5, 0.9), linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = E_N, linetype = "dashed", colour = col_reference) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.25)) +
    labs(
      x = "Total sample size (N)",
      y = "Cumulative probability",
      title = "Cumulative distribution of sample size",
      subtitle = sprintf("Under %s", ss$hypothesis)
    ) +
    theme_publication()

  if (type == "histogram") return(p_hist)
  if (type == "cumulative") return(p_cdf)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_hist + p_cdf + patchwork::plot_layout(ncol = 2)
  } else {
    list(histogram = p_hist, cumulative = p_cdf)
  }
}

#' Plot sample size comparison (H0 vs H1)
plot_sample_size_comparison <- function(ss_compare, type, show_outcomes) {


  dist_h1 <- ss_compare$H1$distribution
  dist_h0 <- ss_compare$H0$distribution

  dist_h1$hypothesis <- "H1"
  dist_h0$hypothesis <- "H0"

  dist_combined <- rbind(dist_h1, dist_h0)
  dist_combined$hypothesis <- factor(dist_combined$hypothesis, levels = c("H1", "H0"))

  # Aggregate
  dist_agg <- aggregate(weight ~ n_total + hypothesis, data = dist_combined, FUN = sum)

  E_N_h1 <- ss_compare$H1$summary$value[ss_compare$H1$summary$metric == "E[N]"]
  E_N_h0 <- ss_compare$H0$summary$value[ss_compare$H0$summary$metric == "E[N]"]

  palette_hyp <- c("H1" = "#2a9d8f", "H0" = "#e76f51")

  bar_width <- max(diff(range(dist_agg$n_total)) / 50, 1)

  p_hist <- ggplot(dist_agg, aes(x = n_total, y = weight, fill = hypothesis)) +
    geom_col(position = "dodge", width = bar_width, alpha = 0.8) +
    geom_vline(xintercept = E_N_h1, linetype = "dashed", colour = palette_hyp["H1"]) +
    geom_vline(xintercept = E_N_h0, linetype = "dashed", colour = palette_hyp["H0"]) +
    scale_fill_manual(values = palette_hyp) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      x = "Total sample size (N)",
      y = "Probability",
      fill = "Hypothesis",
      title = "Sample size distribution",
      subtitle = sprintf("E[N|H1] = %s, E[N|H0] = %s",
                         scales::comma(round(E_N_h1)), scales::comma(round(E_N_h0)))
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  # CDF
  dist_h1_sorted <- dist_h1[order(dist_h1$n_total), ]
  dist_h1_sorted$cum_prob <- cumsum(dist_h1_sorted$weight) / sum(dist_h1_sorted$weight)
  dist_h1_sorted$hypothesis <- "H1"

  dist_h0_sorted <- dist_h0[order(dist_h0$n_total), ]
  dist_h0_sorted$cum_prob <- cumsum(dist_h0_sorted$weight) / sum(dist_h0_sorted$weight)
  dist_h0_sorted$hypothesis <- "H0"

  cdf_combined <- rbind(dist_h1_sorted, dist_h0_sorted)

  p_cdf <- ggplot(cdf_combined, aes(x = n_total, y = cum_prob, colour = hypothesis)) +
    geom_step(linewidth = 0.9) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = palette_hyp) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Total sample size (N)",
      y = "Cumulative probability",
      colour = "Hypothesis",
      title = "Cumulative distribution of sample size"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  if (type == "histogram") return(p_hist)
  if (type == "cumulative") return(p_cdf)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_hist + p_cdf + patchwork::plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    list(histogram = p_hist, cumulative = p_cdf)
  }
}

#' Compare sample sizes under H0 and H1
#' @param results Output from find_lambda_for_power
#' @param sample_size_fn Sample size function
#' @param ... Additional arguments passed to sample_size_fn
compare_sample_sizes <- function(results, sample_size_fn, ...) {

  ss_h1 <- compute_sample_size(results, sample_size_fn, hypothesis = "H1", ...)
  ss_h0 <- compute_sample_size(results, sample_size_fn, hypothesis = "H0", ...)

  comparison <- data.frame(
    metric = ss_h1$summary$metric,
    H1 = ss_h1$summary$value,
    H0 = ss_h0$summary$value
  )

  list(
    comparison = comparison,
    H1 = ss_h1,
    H0 = ss_h0
  )
}

# =============================================================================
# Design Space Plotting
# =============================================================================

#' Plot design space exploration results
#'
#' @param exploration Output from explore_design_space
#' @param type "EN", "Nmax", "clusters", "power", "probabilities", or "all"
#' @param x_var Variable for x-axis (default: "power_stage1")
#' @param group_var Variable for grouping/colour (default: first column of stage1_grid)
plot_design_space <- function(exploration,
                              type = c("EN", "Nmax", "clusters", "power", "probabilities", "all"),
                              x_var = "power_stage1",
                              group_var = NULL) {

  type <- match.arg(type)


  df <- exploration$summary
  stage1_vars <- names(exploration$params$stage1_grid)

  # Default group variable

  if (is.null(group_var)) group_var <- stage1_vars[1]

  df[[group_var]] <- factor(df[[group_var]])
  n_groups <- length(unique(df[[group_var]]))

  palette_groups <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:min(n_groups, 5)]
  if (n_groups > 5) palette_groups <- viridis::viridis(n_groups)

  palette_outcome <- c("Efficacy stop" = "#2d6a4f",
                       "Futility stop" = "#ae2012",
                       "Continue"      = "#457b9d")

  # Label variable (second stage1 var if exists)
  label_var <- if (length(stage1_vars) > 1) stage1_vars[2] else NULL

  # Common plot elements
  make_base_plot <- function(y_var, y_label, title, subtitle = NULL) {
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_line(aes(colour = .data[[group_var]], group = .data[[group_var]]), linewidth = 0.9) +
      geom_point(aes(colour = .data[[group_var]], shape = .data[[group_var]]),
                 size = 3, fill = "white", stroke = 0.8)

    if (!is.null(label_var) && label_var %in% names(df)) {
      p <- p + geom_text(aes(label = .data[[label_var]]),
                         hjust = -0.4, vjust = -0.4, size = 2.8, colour = "grey40")
    }

    p <- p +
      scale_colour_manual(values = palette_groups) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)[1:n_groups]) +
      labs(
        x = if (x_var == "power_stage1") "Stage 1 power" else gsub("_", " ", x_var),
        y = y_label,
        colour = group_var,
        shape = group_var,
        title = title,
        subtitle = subtitle %||% sprintf("Target power = %g%%", exploration$params$target_power * 100)
      ) +
      theme_publication() +
      theme(legend.position = "bottom")

    if (x_var == "power_stage1") {
      p <- p + scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                  expand = expansion(mult = c(0.02, 0.08)))
    }

    p
  }

  # Individual plots
  p_EN <- if ("E_N" %in% names(df)) {
    make_base_plot("E_N", expression(E*"[N]"), "Expected total sample size") +
      scale_y_continuous(labels = scales::comma)
  } else NULL

  p_Nmax <- if ("N_max" %in% names(df)) {
    make_base_plot("N_max", expression(N[max]), "Maximum total sample size",
                   "Worst-case if trial continues to stage 2") +
      scale_y_continuous(labels = scales::comma)
  } else NULL

  p_clusters <- if ("E_clusters" %in% names(df)) {
    make_base_plot("E_clusters", "E[total clusters]", "Expected total clusters",
                   "Including new clusters in stage 2")
  } else if ("E_total_clusters" %in% names(df)) {
    make_base_plot("E_total_clusters", "E[total clusters]", "Expected total clusters")
  } else NULL

  p_power <- if (length(stage1_vars) > 1 && "power_stage1" %in% names(df)) {
    x2 <- stage1_vars[2]
    ggplot(df, aes(x = .data[[x2]], y = power_stage1)) +
      geom_line(aes(colour = .data[[group_var]], group = .data[[group_var]]), linewidth = 0.9) +
      geom_point(aes(colour = .data[[group_var]], shape = .data[[group_var]]),
                 size = 3, fill = "white", stroke = 0.8) +
      scale_colour_manual(values = palette_groups) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)[1:n_groups]) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = gsub("_", " ", x2),
        y = "Stage 1 power",
        colour = group_var,
        shape = group_var,
        title = "Stage 1 power by design parameters"
      ) +
      theme_publication() +
      theme(legend.position = "bottom")
  } else NULL

  # Probabilities plot
  p_probs <- if (all(c("prob_efficacy", "prob_futility", "prob_continue") %in% names(df))) {
    df_long <- tidyr::pivot_longer(
      df,
      cols = c(prob_efficacy, prob_futility, prob_continue),
      names_to = "outcome",
      values_to = "probability",
      names_prefix = "prob_"
    )
    df_long$outcome <- factor(
      df_long$outcome,
      levels = c("efficacy", "continue", "futility"),
      labels = c("Efficacy stop", "Continue", "Futility stop")
    )

    ggplot(df_long, aes(x = .data[[x_var]], y = probability, colour = outcome)) +
      geom_line(aes(linetype = outcome), linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(as.formula(paste("~", group_var)),
                 labeller = labeller(.default = function(x) paste(group_var, "=", x))) +
      scale_colour_manual(values = palette_outcome) +
      scale_linetype_manual(values = c("Efficacy stop" = "solid",
                                       "Continue" = "dashed",
                                       "Futility stop" = "dotted")) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
      labs(
        x = if (x_var == "power_stage1") "Stage 1 power" else gsub("_", " ", x_var),
        y = "Probability",
        colour = "Outcome",
        linetype = "Outcome",
        title = "Interim decision probabilities"
      ) +
      theme_publication() +
      theme(strip.text = element_text(face = "bold"),
            legend.position = "bottom")
  } else NULL

  # Return requested plot(s)
  if (type == "EN") return(p_EN)
  if (type == "Nmax") return(p_Nmax)
  if (type == "clusters") return(p_clusters)
  if (type == "power") return(p_power)
  if (type == "probabilities") return(p_probs)

  if (type == "all") {
    plots <- Filter(Negate(is.null), list(EN = p_EN, Nmax = p_Nmax,
                                          clusters = p_clusters, power = p_power))
    if (requireNamespace("patchwork", quietly = TRUE) && length(plots) > 1) {
      patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    } else {
      plots
    }
  }
}

#' Create heatmap of design space
#'
#' @param exploration Output from explore_design_space
#' @param z_var Variable for fill colour
#' @param x_var, y_var Variables for axes (default: first two stage1 vars)
plot_design_heatmap <- function(exploration,
                                z_var = "E_N",
                                x_var = NULL,
                                y_var = NULL) {


  df <- exploration$summary
  stage1_vars <- names(exploration$params$stage1_grid)

  if (is.null(x_var)) x_var <- stage1_vars[1]
  if (is.null(y_var) && length(stage1_vars) > 1) y_var <- stage1_vars[2]

  if (is.null(y_var)) {
    stop("Need at least 2 stage 1 variables for heatmap")
  }

  z_labels <- c(
    E_N = "E[N]",
    N_max = "Maximum N",
    E_total_clusters = "E[clusters]",
    E_clusters = "E[clusters]",
    lambda = "Optimal λ",
    power_stage1 = "Stage 1 power"
  )

  z_label <- z_labels[z_var] %||% z_var

  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_tile(aes(fill = .data[[z_var]])) +
    geom_text(aes(label = sprintf("%.0f", .data[[z_var]])), colour = "white", size = 3) +
    scale_fill_viridis_c(option = "plasma") +
    labs(x = gsub("_", " ", x_var),
         y = gsub("_", " ", y_var),
         fill = z_label,
         title = sprintf("%s across design space", z_label)) +
    theme_publication()
}

# =============================================================================
# Fixed Design Comparison
# =============================================================================

#' Find optimal fixed (non-adaptive) design
#'
#' @param target_power Target power
#' @param design_grid Grid of fixed design parameters (e.g., expand.grid(k=..., m=...))
#' @param model_fn_factory Factory for model function
#' @param fixed_params_base Fixed parameters
#' @param power_fn Function to compute power from model output.
#'        Signature: function(model_output) -> power (scalar)
#'        Default uses I_eff and b1 for normal approximation.
#' @param cost_fn Cost function for fixed design
#' @param cost_params Cost parameters
#' @param verbose Print progress
#' @return List with optimal design and all evaluated designs
find_optimal_fixed_design <- function(target_power = 0.8,
                                      design_grid,
                                      model_fn_factory,
                                      fixed_params_base,
                                      power_fn = NULL,
                                      cost_fn,
                                      cost_params = list(rho = 25),
                                      verbose = TRUE) {

  # Default power function: two-sided test using I_eff
  if (is.null(power_fn)) {
    power_fn <- function(model_output) {
      mu <- model_output$b1 * sqrt(model_output$I_eff)
      pnorm(-1.96 - mu) + (1 - pnorm(1.96 - mu))
    }
  }

  n_designs <- nrow(design_grid)
  if (verbose) cat(sprintf("Evaluating %d fixed designs...\n", n_designs))

  results <- lapply(1:n_designs, function(i) {
    design_params <- design_grid[i, , drop = FALSE]

    # Merge with fixed params
    fixed_params <- fixed_params_base
    for (nm in names(design_params)) {
      fixed_params[[nm]] <- design_params[[nm]]
    }

    # Get model output
    model_fn <- model_fn_factory()

    # For fixed design, stage 2 params = stage 1 params (single stage)
    model_output <- tryCatch({
      model_fn(design_params, fixed_params)
    }, error = function(e) NULL)

    if (is.null(model_output)) {
      return(data.frame(design_params, power = NA, cost = NA))
    }

    power <- power_fn(model_output)
    cost <- do.call(cost_fn, c(list(model_output), cost_params))

    data.frame(design_params, power = power, cost = cost)
  })

  results_df <- do.call(rbind, results)

  # Find designs meeting power target
  valid <- results_df$power >= target_power & !is.na(results_df$power)

  if (!any(valid)) {
    warning("No designs achieve target power")
    optimal <- results_df[which.max(results_df$power), ]
  } else {
    # Among valid, find minimum cost
    valid_df <- results_df[valid, ]
    optimal <- valid_df[which.min(valid_df$cost), ]
  }

  if (verbose) {
    cat(sprintf("Optimal fixed design: power = %.3f, cost = %.1f\n",
                optimal$power, optimal$cost))
  }

  list(
    optimal = optimal,
    all_designs = results_df,
    target_power = target_power
  )
}

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
      # Generic fallback
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
      rho = rho
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
    # Parallel: vary k and m
    k_vals <- unique(stage1_grid$k1)
    m_vals <- unique(stage1_grid$m1)

    grid <- expand.grid(
      k = seq(min(k_vals), ceiling(max(k_vals) * 1.5), by = 1),
      m = seq(min(m_vals), ceiling(max(m_vals) * 2), by = 5)
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

    grid <- expand.grid(
      k = seq(min(k_vals), ceiling(max(k_vals) * 1.5), by = 1),
      m = seq(min(m_vals), ceiling(max(m_vals) * 2), by = 5)
    )

    model_fn <- fixed_crossover_model_builder

    cost_fn <- function(model_output, rho) {
      n_cost <- model_output$n_total
      k_cost <- model_output$k * 2
      n_cost + rho * k_cost
    }

  } else if (spec$design_type == "stepped_wedge") {
    # Stepped-wedge: K is typically fixed, vary T (total periods) and possibly m
    K <- unique(stage1_grid$K)[1]
    t1_vals <- unique(stage1_grid$t1)
    m_val <- unique(stage1_grid$m)[1]

    # For fixed SW, we vary total time periods
    # The "fixed" comparison is a complete SW with no interim
    grid <- expand.grid(
      K = K,
      T_total = seq(max(t1_vals), max(t1_vals) * 3, by = 1),
      m = m_val
    )

    model_fn <- fixed_stepped_wedge_model_builder

    cost_fn <- function(model_output, rho) {
      n_cost <- model_output$n_total
      t_cost <- model_output$T_total
      n_cost + rho * t_cost
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

#' Fixed stepped-wedge model builder
#' @keywords internal
fixed_stepped_wedge_model_builder <- function(design_params, fixed_params) {

  K <- design_params$K
  T_total <- design_params$T_total
  m <- design_params$m
  icc <- fixed_params$icc
  cac <- fixed_params$cac %||% 0.8
  delta <- fixed_params$delta


  # Standard SW: one cluster switches per period (or spread evenly)
  r <- max(1, floor(K / (T_total - 1)))  # Rollout rate

  # Calculate switch times
  switch_times <- numeric(K)
  clusters_assigned <- 0
  for (t in 2:T_total) {  # First period is baseline (all control)
    n_switch <- min(r, K - clusters_assigned)
    if (n_switch > 0) {
      idx <- (clusters_assigned + 1):(clusters_assigned + n_switch)
      switch_times[idx] <- t
      clusters_assigned <- clusters_assigned + n_switch
    }
    if (clusters_assigned >= K) break
  }
  # Any remaining switch at last period
  if (any(switch_times == 0)) {
    switch_times[switch_times == 0] <- T_total
  }

  # Build data frame
  df <- expand.grid(cl = 1:K, t = 1:T_total)
  df$trt <- as.integer(df$t >= switch_times[df$cl])
  df$n <- m

  v1 <- c(icc * cac, icc * (1 - cac))

  mod <- Model$new(
    as.formula("~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t))"),
    data = df,
    family = gaussian(),
    mean = c(rep(0, T_total), delta, 0),
    covariance = v1,
    weights = df$n,
    var_par = 1 - icc
  )

  X <- mod$mean$X
  D <- mod$covariance$D
  Z <- mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  V_inv_X <- solve(V, X)
  I_full <- crossprod(X, V_inv_X)

  # Treatment is the last fixed effect before random
  j_trt <- ncol(X) - 1
  I_eff <- 1 / solve(I_full)[j_trt, j_trt]

  n_total <- K * T_total * m

  list(
    I_eff = I_eff,
    b1 = delta,
    K = K,
    T_total = T_total,
    m = m,
    r = r,
    n_total = n_total
  )
}

#' Fixed crossover model builder
#' @keywords internal
fixed_crossover_model_builder <- function(design_params, fixed_params) {

  k <- design_params$k
  m <- design_params$m
  icc <- fixed_params$icc
  cac <- fixed_params$cac %||% 0.8
  delta <- fixed_params$delta

  # Two-period crossover
  df <- nelder(as.formula(paste0("~ int(2) > cl(", k, ") > t(2)")))
  df$int <- df$int - 1
  df$t <- df$t - 1
  df$trt <- ifelse((df$int == 0 & df$t == 1) | (df$int == 1 & df$t == 0), 1, 0)
  df$n <- m

  v1 <- c(icc * cac, icc * (1 - cac))

  mod <- Model$new(
    ~ factor(t) + trt + (1|gr(cl)) + (1|gr(cl,t)),
    data = df,
    family = gaussian(),
    mean = c(0, 0, delta, 0),
    covariance = v1,
    weights = df$n,
    var_par = 1 - icc
  )

  X <- mod$mean$X
  D <- mod$covariance$D
  Z <- mod$covariance$Z

  Z_sp <- Matrix(Z, sparse = TRUE)
  D_sp <- Matrix(D, sparse = TRUE)
  S <- Diagonal(x = 1/mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  V_inv_X <- solve(V, X)
  I_full <- crossprod(X, V_inv_X)
  j_trt <- ncol(X) - 1
  I_eff <- 1 / solve(I_full)[j_trt, j_trt]

  list(
    I_eff = I_eff,
    b1 = delta,
    k = k,
    m = m,
    n_total = 2 * k * m * 2  # 2 arms, 2 periods
  )
}

#' Plot adaptive vs fixed design comparison
#'
#' @param comparison Output from compare_to_fixed
#' @param type "cost_ratio", "savings", "frontier", or "all"
plot_fixed_comparison <- function(comparison,
                                  type = c("cost_ratio", "savings", "frontier", "all")) {

  type <- match.arg(type)


  df <- comparison$comparison
  fixed_cost <- comparison$fixed$optimal$cost

  stage1_vars <- names(comparison$adaptive$params$stage1_grid)
  group_var <- stage1_vars[1]
  df[[group_var]] <- factor(df[[group_var]])
  n_groups <- length(unique(df[[group_var]]))

  palette_groups <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[1:n_groups]

  # Cost ratio plot
  p_ratio <- ggplot(df, aes(x = power_stage1)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(aes(y = E_cost_ratio, colour = .data[[group_var]]), linewidth = 0.9) +
    geom_line(aes(y = max_cost_ratio, colour = .data[[group_var]]),
              linetype = "dotted", linewidth = 0.7) +
    geom_point(aes(y = E_cost_ratio, colour = .data[[group_var]]), size = 2) +
    scale_colour_manual(values = palette_groups) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Stage 1 power",
      y = "Cost ratio (adaptive / fixed)",
      colour = group_var,
      title = "Cost comparison: adaptive vs fixed",
      subtitle = "Solid = E[cost], Dotted = max cost"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  # Savings plot
  p_savings <- ggplot(df, aes(x = power_stage1, y = E_savings_pct)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(aes(colour = .data[[group_var]]), linewidth = 0.9) +
    geom_point(aes(colour = .data[[group_var]], shape = max_within_fixed), size = 3) +
    scale_colour_manual(values = palette_groups) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                       labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       name = "Max ≤ fixed") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Stage 1 power",
      y = "Expected savings (%)",
      colour = group_var,
      title = "Expected cost savings vs fixed design"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  # Frontier plot
  p_frontier <- ggplot(df, aes(x = E_cost, y = max_cost)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = fixed_cost, linetype = "dashed", colour = "#e76f51") +
    geom_hline(yintercept = fixed_cost, linetype = "dashed", colour = "#e76f51") +
    geom_point(aes(colour = .data[[group_var]], shape = max_within_fixed), size = 3) +
    scale_colour_manual(values = palette_groups) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
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
# Pareto Summary
# =============================================================================

#' Summarise Pareto frontier
#'
#' @param pareto_df Output from find_constrained_pareto
#' @return Summary data frame (invisibly)
summarise_pareto <- function(pareto_df) {

  obj_list <- attr(pareto_df, "objectives")

  if (is.null(obj_list)) {
    message("No objectives attribute found; returning basic summary")
    return(summary(pareto_df))
  }

  obj_names <- names(obj_list)
  obj_directions <- unlist(obj_list)

  summary_list <- lapply(seq_along(obj_names), function(i) {
    col <- obj_names[i]
    dir <- obj_directions[i]
    vals <- pareto_df[[col]]

    data.frame(
      objective = col,
      direction = dir,
      min = min(vals, na.rm = TRUE),
      max = max(vals, na.rm = TRUE),
      range = max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE),
      best = if (dir == "min") min(vals, na.rm = TRUE) else max(vals, na.rm = TRUE)
    )
  })

  summary_df <- do.call(rbind, summary_list)

  cat(sprintf("Pareto frontier: %d designs\n", nrow(pareto_df)))
  cat(sprintf("Dominated designs excluded: %d\n",
              attr(pareto_df, "n_dominated") %||% NA))
  cat("\nObjective ranges on frontier:\n")
  print(summary_df, row.names = FALSE)

  invisible(summary_df)
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
