sw_model_builder <- function(design_params, fixed_params) {


  # Access cache from fixed_params (ISSUE 1: was undefined)
  cache <- fixed_params$cache

  generate_sw_design <- function(k, t1, t2, r) {
    cl_mid <- (k - t1 + 1) / 2

    # Stage 1
    df1 <- glmmrBase::nelder(as.formula(paste0("~ cl(", k, ") * t(", t1, ")")))
    df1$int <- I(df1$t > df1$cl) * 1

    # Stage 2
    df2 <- glmmrBase::nelder(as.formula(paste0("~ cl(", (k - t1 + 1), ") * t(", t2, ")")))
    df3 <- glmmrBase::nelder(as.formula(paste0("~ cl(", t1 - 1, ") * t(", t2, ")")))
    df3$int <- 1
    df3$t <- df3$t + t1

    # Intervention indicator for stage 2 middle clusters
    if (r < 1e-6) {
      df2$int <- (df2$cl <= cl_mid) * 1
    } else {
      df2$int <- (df2$t > cl_mid + (df2$cl - cl_mid) / r) * 1
    }
    df2$t <- df2$t + t1
    df2$cl <- df2$cl + (t1 - 1)

    # Offset df3 cluster IDs to avoid collision
    # df3$cl <- df3$cl + k  # ISSUE 2: df3 clusters need offset

    df <- rbind(df1, df2, df3)
    return(df)
  }

  # Extract parameters
  k <- fixed_params$k
  t1 <- fixed_params$t1
  t2 <- design_params$t2
  r <- design_params$r
  m1 <- fixed_params$m1
  m2 <- design_params$m2

  # ISSUE 3: Function call had wrong number of arguments
  df <- generate_sw_design(k, t1, t2, r)

  # Set weights
  df$n <- m1
  df[df$t > t1, 'n'] <- m2  # ISSUE 4: t1 was undefined, now using extracted variable

  # Convert ICC to variance parameter for binomial
  icc_var_par <- function(baseline, icc) {
    v1 <- baseline * (1 - baseline)
    v2 <- v1 * (icc / (1 - icc))
    tausq <- v2 / (baseline^2 / ((1 + exp(log(baseline / (1 - baseline))))^2))
    return(tausq)
  }

  to_logit_scale <- function(baseline, eff_size) {
    b0 <- log(baseline / (1 - baseline))
    m1 <- baseline + eff_size
    b1 <- log(m1 / (1 - m1)) - b0
    return(c(b0, b1))
  }

  beta <- to_logit_scale(fixed_params$baseline, fixed_params$delta)
  v1 <- icc_var_par(fixed_params$baseline, fixed_params$icc)

  # ISSUE 5: ar was undefined - need to define or get from fixed_params
  ar <- fixed_params$cac  # Assuming cac is the AR correlation parameter

  T_total <- t1 + t2

  # ISSUE 6: cache_key used undefined variables
  cache_key <- paste(k, t1, t2, r, sep = "_")

  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- glmmrBase::Model$new(
      ~ int + factor(t) + (1|gr(cl)*ar0(t)),
      data = df,
      family = binomial(),
      mean = c(beta, rep(0, T_total - 1)),  # intercept handled by beta[1], then trt, then period dummies
      covariance = c(v1, ar),
      trials = df$n
    )
    cache$cache_key <- cache_key
  } else {
    glmmrBase:::Model__set_trials(cache$mod$.__enclos_env__$private$ptr, df$n)
    cache$mod$update_parameters(cov.pars = c(v1, ar))
  }

  n1 <- nrow(df[df$t <= t1, ])  # More robust than assuming df1 exists here
  n2 <- nrow(df) - n1

  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z

  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)

  idx1 <- which(df$t <= t1)  # More robust indexing
  idx2 <- which(df$t > t1)

  # ISSUE 7: j_trt calculation - int is column 2 in ~ int + factor(t) + ...
  # Column 1 = intercept (implicit), Column 2 = int (treatment)
  j_trt <- 2

  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)

  df_s1 <- t1 * k - t1 - 1
  df_full <- (t1 + t2) * k - t1 - t2 - 1

  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = beta[2],  # ISSUE 8: Should be treatment effect (beta[2]), not intercept (beta[1])
    n1 = n1,
    n2 = n2,
    m1 = m1,
    m2 = m2,
    N1_total = k * t1 * m1,       # individuals in stage 1
    N2_total = k * t2 * m2,       # individuals in stage 2
    N_total = k * t1 * m1 + k * t2 * m2,  # total individuals
    k = k,
    t1 = t1,
    t2 = t2,
    r = r,
    T_total = T_total,
    df_s1 = df_s1,
    df_full = df_full
  )
}

# Also need factory function for caching
make_sw_model_fn <- function() {
  cache <- new.env()
  cache$mod <- NULL
  cache$cache_key <- NULL

  function(design_params, fixed_params) {
    fixed_params$cache <- cache
    sw_model_builder(design_params, fixed_params)
  }
}

sw_crt <- function(icc, cac = 0.8, delta, baseline,
                   k, m1, t1,
                   t2 = 1:10,
                   r = seq(0, 1.2, by = 0.2),
                   m2 = seq(20, 100, by = 20),
                   rho = 30) {

  spec <- crt_design_spec(
    stage1_params = c("k", "t1", "m1"),
    stage2_params = c("t2", "r", "m2"),
    resources = list(
      n_s1 = ~ k * t1 * m1,
      n_s2 = ~ k * t2 * m2,
      periods_s1 = ~ k * t1,
      periods_s2 = ~ k * t2
    ),
    cost_structure = list(
      weights = c(n = 1, periods = "rho"),
      stage2_resources = c("n_s2", "periods_s2")
    ),
    model_builder = sw_model_builder,
    n_arms = 1,
    design_type = "stepped_wedge"
  )

  structure(
    list(
      spec = spec,
      fixed_params = list(
        icc = icc,
        cac = cac,
        delta = delta,
        baseline = baseline,
        k = k
      ),
      stage1_grid = expand.grid(t1 = t1, m1 = m1),
      stage2_grid_fn = function(s1) expand.grid(t2 = t2, r = r, m2 = m2),
      rho = rho
    ),
    class = c("sw_crt", "crt_design")
  )
}


# Define the design in ONE call - no custom functions needed!
design <- sw_crt(
  icc = 0.06,           # Intra-cluster correlation
  cac = 0.8,            # Cluster autocorrelation
  delta = -0.07,         # Treatment effect
  baseline = 0.2,
  k = 11,    # Stage 1 clusters per arm (explore these)
  t1 = 3:6,
  m1 = seq(10,60,by=10),
  rho = 15              # Cluster-to-individual cost ratio
)

# Run the analysis - ONE function call
results <- adaptive_analysis(design, target_power = 0.8, verbose = TRUE, tol = 0.01, method = "cost_cap")

# View summary
summary(results)

# Plot results
p1 <- plot(results, type = "EN")      # Expected sample size
p2 <- plot(results, type = "Nmax")    # Maximum sample size
p3 <- plot(results, type = "pareto", objectives = list(E_cost = "min", max_cost = "min"))  # Pareto frontier

pareto <- find_pareto(results, objectives = list(E_cost = "min", max_cost = "min"))


# Define the design in ONE call - no custom functions needed!
design_single <- sw_crt(
  icc = 0.06,           # Intra-cluster correlation
  cac = 0.8,            # Cluster autocorrelation
  delta = -0.07,         # Treatment effect
  baseline = 0.2,
  k = 11,    # Stage 1 clusters per arm (explore these)
  t1 = 5, # Stage 1 individuals per cluster
  m1 = 30,
  rho = 15              # Cluster-to-individual cost ratio
)

# Analyse
results_single <- adaptive_analysis(design_single, target_power = 0.8, tol = 0.01, method = "cost_cap")

interim_analysis(results_single, z1_obs = -1.7,
                 theta_hat = list(icc = 0.001))

int1 <-interim_sensitivity(
  results_single,
  scenarios = list(
    "Planned (ICC = 0.06)" = list(icc = 0.06),
    "Favourable (ICC = 0.02)"           = list(icc = 0.02),
    "Very favourable (ICC = 0.01)"         = list(icc = 0.01)
  )
)

# Summary shows more detail for single designs
summary(int1)

plot(int1, type = "decision")
plot(int1, type = "cp")
plot(sens, type = "all")

# Decision rules
plot(results_single, type = "decision")
plot_decision_rules(results_single, design_var = "t2")
# Get decision rules as data frame
rules <- get_decision_rules(results_single)
print(rules[seq(1, nrow(rules), by = 5), c("z1", "continue", "cp", "m2", "k2")])

plot_sw_design <- function(df, t1 = NULL, title = NULL,
                           show_legend = TRUE,
                           cluster_labels = NULL,
                           time_labels = NULL) {


  # Ensure factors for proper ordering
  df$cl <- factor(df$cl, levels = sort(unique(df$cl), decreasing = TRUE))
  df$t <- factor(df$t, levels = sort(unique(df$t)))
  df$int <- factor(df$int, levels = c(0, 1), labels = c("Control", "Intervention"))

  # Colour palette
  palette_trt <- c("Control" = "#f0f0f0", "Intervention" = "#2a9d8f")

  # Base plot
  p <- ggplot(df, aes(x = t, y = cl, fill = int)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = palette_trt, name = "Status") +
    labs(
      x = "Time period",
      y = "Cluster",
      title = title %||% "Stepped-Wedge Trial Design"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )

  # Add interim analysis line
  if (!is.null(t1)) {
    # Position line between t1 and t1+1
    p <- p +
      geom_vline(xintercept = t1 + 0.5, colour = "#e76f51",
                 linewidth = 1.2, linetype = "dashed") +
      annotate("text", x = t1 + 0.5, y = length(unique(df$cl)) + 0.7,
               label = "Interim", colour = "#e76f51",
               fontface = "bold", size = 3.5, hjust = 0.5)
  }

  # Custom axis labels if provided
  if (!is.null(cluster_labels)) {
    p <- p + scale_y_discrete(labels = rev(cluster_labels))
  }
  if (!is.null(time_labels)) {
    p <- p + scale_x_discrete(labels = time_labels)
  }

  p
}

#' Plot stepped-wedge design with sample size weights
#'
#' @param df Data frame with columns: cl, t, int, and optionally n (sample size)
#' @param t1 Stage 1 periods
#' @param show_n Show sample sizes in cells (default TRUE)
#' @param title Optional title
#' @return A ggplot object
plot_sw_design_detailed <- function(df, t1 = NULL, show_n = TRUE, title = NULL) {


  # Ensure factors
  df$cl <- factor(df$cl, levels = sort(unique(df$cl), decreasing = TRUE))
  df$t <- factor(df$t, levels = sort(unique(df$t)))
  df$int_label <- factor(df$int, levels = c(0, 1), labels = c("Control", "Intervention"))

  # Colour palette
  palette_trt <- c("Control" = "#e8e8e8", "Intervention" = "#2a9d8f")

  p <- ggplot(df, aes(x = t, y = cl)) +
    geom_tile(aes(fill = int_label), colour = "white", linewidth = 0.8) +
    scale_fill_manual(values = palette_trt, name = "Status")

  # Add sample sizes if present and requested
  if (show_n && "n" %in% names(df)) {
    p <- p + geom_text(aes(label = n), colour = "grey30", size = 3)
  }

  p <- p +
    labs(
      x = "Time period",
      y = "Cluster",
      title = title %||% "Stepped-Wedge Trial Design"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )

  # Interim analysis line
  if (!is.null(t1)) {
    n_clusters <- length(unique(df$cl))
    p <- p +
      geom_vline(xintercept = t1 + 0.5, colour = "#e76f51",
                 linewidth = 1.5, linetype = "dashed") +
      annotate("label", x = t1 + 0.5, y = n_clusters + 0.8,
               label = "Interim Analysis", fill = "#e76f51", colour = "white",
               fontface = "bold", size = 3, label.padding = unit(0.2, "lines"))
  }

  # Add stage labels
  if (!is.null(t1)) {
    t_max <- max(as.numeric(as.character(df$t)))
    p <- p +
      annotate("text", x = (1 + t1) / 2, y = 0.3,
               label = "Stage 1", fontface = "italic", size = 3.5, colour = "grey40") +
      annotate("text", x = t1 + (t_max - t1) / 2 + 0.5, y = 0.3,
               label = "Stage 2", fontface = "italic", size = 3.5, colour = "grey40")
  }

  p + coord_cartesian(clip = "off", ylim = c(0.5, length(unique(df$cl)) + 0.5))
}

#' Plot multiple SW designs side by side for comparison
#'
#' @param designs Named list of data frames
#' @param t1 Stage 1 periods (same for all)
#' @return A combined ggplot
plot_sw_designs_compare <- function(designs, t1 = NULL) {

  # Combine with design identifier
  df_list <- lapply(names(designs), function(nm) {
    d <- designs[[nm]]
    d$design <- nm
    d
  })
  df_all <- do.call(rbind, df_list)

  df_all$cl <- factor(df_all$cl, levels = sort(unique(df_all$cl), decreasing = TRUE))
  df_all$t <- factor(df_all$t, levels = sort(unique(df_all$t)))
  df_all$int_label <- factor(df_all$int, levels = c(0, 1), labels = c("Control", "Intervention"))
  df_all$design <- factor(df_all$design, levels = names(designs))

  palette_trt <- c("Control" = "#e8e8e8", "Intervention" = "#2a9d8f")

  p <- ggplot(df_all, aes(x = t, y = cl, fill = int_label)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = palette_trt, name = "Status") +
    facet_wrap(~ design, ncol = length(designs)) +
    labs(x = "Time period", y = "Cluster") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )

  if (!is.null(t1)) {
    p <- p + geom_vline(xintercept = t1 + 0.5, colour = "#e76f51",
                        linewidth = 1, linetype = "dashed")
  }

  p
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

generate_sw_design <- function(k, t1, t2, r) {
  cl_mid <- (k - t1 + 1) / 2

  # Stage 1
  df1 <- glmmrBase::nelder(as.formula(paste0("~ cl(", k, ") * t(", t1, ")")))
  df1$int <- I(df1$t > df1$cl) * 1

  # Stage 2
  df2 <- glmmrBase::nelder(as.formula(paste0("~ cl(", (k - t1 + 1), ") * t(", t2, ")")))
  df3 <- glmmrBase::nelder(as.formula(paste0("~ cl(", t1 - 1, ") * t(", t2, ")")))
  df3$int <- 1
  df3$t <- df3$t + t1

  # Intervention indicator for stage 2 middle clusters
  if (r < 1e-6) {
    df2$int <- (df2$cl <= cl_mid) * 1
  } else {
    df2$int <- (df2$t > cl_mid + (df2$cl - cl_mid) / r) * 1
  }
  df2$t <- df2$t + t1
  df2$cl <- df2$cl + (t1 - 1)

  # Offset df3 cluster IDs to avoid collision
  # df3$cl <- df3$cl + k  # ISSUE 2: df3 clusters need offset

  df <- rbind(df1, df2, df3)
  return(df)
}

df <- generate_sw_design(11,4,5,1)
plot_sw_design(df, t1 = 4, title = "SW Design: k=11, t1=4, t2=5, r=1.0")


df$n <- ifelse(df$t <= 4, 20, 50)

designs <- list(
  "ICC = 0.06 (m2 = 80)" = generate_sw_design(11, 5, 7, r = 1.05),
  "ICC = 0.02 (m2 = 40)" = generate_sw_design(11, 5, 3, r = 0.7),
  "ICC = 0.01 (m2 = 60)" = generate_sw_design(11, 5, 1, r = 0.0)
)
plot_sw_designs_compare(designs, t1 = 4)

