# =============================================================================
# Design Specification System
# =============================================================================
# Provides a declarative way to specify CRT designs, from which cost functions,
# sample size functions, and other components are auto-generated.

#' Create a CRT design specification
#'
#' @param stage1_params Character vector of stage 1 parameter names
#' @param stage2_params Character vector of stage 2 parameter names
#' @param resources Named list of formulas for resource calculation.
#'        Names should follow convention: {resource_type}_{stage}
#'        e.g., n_s1, n_s2, clusters_s1, clusters_s2, periods_s2
#' @param cost_structure Named list defining cost components:
#'        - weights: named vector of cost weights (can reference 'rho')
#'        - stage2_only: which resources are stage 2 only (for λ-calibration)
#' @param model_builder Function to build the statistical model
#' @param n_arms Number of arms (default 2)
#' @param design_type Character string identifying the design type
#' @return A crt_design_spec object
crt_design_spec <- function(stage1_params,
                            stage2_params,
                            resources,
                            cost_structure = NULL,
                            model_builder,
                            n_arms = 2,
                            design_type = "custom") {


  # Default cost structure
  if (is.null(cost_structure)) {
    cost_structure <- list(
      weights = c(n = 1, clusters = "rho", periods = "rho"),
      stage2_resources = grep("_s2$|_stage2$|2$", names(resources), value = TRUE)
    )
  }

  structure(
    list(
      stage1_params = stage1_params,
      stage2_params = stage2_params,
      resources = resources,
      cost_structure = cost_structure,
      model_builder = model_builder,
      n_arms = n_arms,
      design_type = design_type
    ),
    class = "crt_design_spec"
  )
}

#' Evaluate a resource formula in an environment
#' @param formula A formula like ~ n_arms * k1 * m1
#' @param env A list or environment with variable values
#' @return Numeric result
eval_resource <- function(formula, env) {
  if (is.null(formula)) return(0)
  if (inherits(formula, "formula")) {
    eval(formula[[2]], envir = as.list(env))
  } else if (is.numeric(formula)) {
    formula
  } else {
    eval(parse(text = as.character(formula)), envir = as.list(env))
  }
}

#' Get cost weight, resolving references to variables
#' @param weight_spec Weight specification (number or variable name as string)
#' @param env Environment with variables
#' @return Numeric weight
resolve_weight <- function(weight_spec, env) {
  if (is.numeric(weight_spec)) {
    weight_spec
  } else if (is.character(weight_spec)) {
    env[[weight_spec]] %||% 1

} else {
    1
  }
}

#' Identify which resource type a resource name belongs to
#' @param resource_name Name like "n_s1", "clusters_s2"
#' @param cost_weights Named vector of weights like c(n = 1, clusters = "rho")
#' @return The matching weight name or NULL
match_resource_type <- function(resource_name, cost_weights) {
  for (type in names(cost_weights)) {
    if (grepl(paste0("^", type), resource_name, ignore.case = TRUE)) {
      return(type)
    }
  }
  # Try without prefix
  base_name <- sub("_s[12]$|_stage[12]$", "", resource_name)
  if (base_name %in% names(cost_weights)) {
    return(base_name)
  }
  NULL
}

# =============================================================================
# Cost Function Generation
# =============================================================================

#' Generate cost function from specification
#'
#' @param spec A crt_design_spec object
#' @param stage Which stage cost: "stage2" (for λ-calibration) or "total"
#' @return A function with signature function(model_output, rho)
generate_cost_fn <- function(spec, stage = c("stage2", "total")) {
  stage <- match.arg(stage)

  # Capture spec in closure
  resources <- spec$resources
  cost_structure <- spec$cost_structure
  n_arms <- spec$n_arms

  function(model_output, rho) {
    # Build evaluation environment
    env <- if (is.list(model_output)) as.list(model_output) else list()
    env$rho <- rho
    env$n_arms <- n_arms

    # Determine which resources to include
    if (stage == "stage2") {
      resource_names <- cost_structure$stage2_resources
      if (length(resource_names) == 0) {
        # Fallback: resources with s2/stage2/2 in name
        resource_names <- grep("s2|stage2|_2$", names(resources),
                               value = TRUE, ignore.case = TRUE)
      }
    } else {
      resource_names <- names(resources)
    }

    # Calculate total cost
    total_cost <- 0
    for (rname in resource_names) {
      if (!rname %in% names(resources)) next

      # Evaluate resource quantity
      qty <- tryCatch(
        eval_resource(resources[[rname]], env),
        error = function(e) 0
      )

      # Find matching weight
      weight_type <- match_resource_type(rname, cost_structure$weights)
      weight <- if (!is.null(weight_type)) {
        resolve_weight(cost_structure$weights[[weight_type]], env)
      } else {
        1
      }

      total_cost <- total_cost + qty * weight
    }

    total_cost
  }
}

# =============================================================================
# Sample Size Function Generation
# =============================================================================

#' Generate sample size function from specification
#'
#' @param spec A crt_design_spec object
#' @return A function with signature function(opt_designs, results, rho)
generate_sample_size_fn <- function(spec) {

  # Capture spec components
  resources <- spec$resources
  cost_structure <- spec$cost_structure
  stage2_params <- spec$stage2_params
  n_arms <- spec$n_arms

  # Find the n resources
  n_s1_name <- grep("^n_s1$|^n_stage1$", names(resources), value = TRUE)[1]
  n_s2_name <- grep("^n_s2$|^n_stage2$", names(resources), value = TRUE)[1]

  # Find cluster resources (if any)
  clusters_s1_name <- grep("^clusters_s1|^k_s1", names(resources), value = TRUE)[1]
  clusters_s2_name <- grep("^clusters_s2|^k_s2", names(resources), value = TRUE)[1]

  function(opt_designs, results, rho = 30) {

    # Get model output for base parameters
    model1 <- results$models$list[[1]]
    design_grid <- results$models$design_grid

    # Base environment with stage 1 params
    env_base <- as.list(model1)
    env_base$rho <- rho
    env_base$n_arms <- n_arms

    # Stage 1 sample size (fixed)
    n_stage1 <- if (!is.na(n_s1_name)) {
      eval_resource(resources[[n_s1_name]], env_base)
    } else {
      # Fallback
      n_arms * (env_base$k1 %||% 1) * (env_base$m1 %||% 1)
    }

    # Stage 2 sample sizes (vector, one per z1 point)
    n_stage2 <- sapply(seq_len(nrow(opt_designs)), function(i) {
      if (!opt_designs$continue[i]) return(0)

      env <- env_base
      for (p in stage2_params) {
        if (p %in% names(opt_designs)) {
          val <- opt_designs[[p]][i]
          env[[p]] <- if (is.na(val)) 0 else val
        }
      }

      if (!is.na(n_s2_name)) {
        eval_resource(resources[[n_s2_name]], env)
      } else {
        0
      }
    })

    # Maximum stage 2 (from design grid maxima)
    env_max <- env_base
    for (p in stage2_params) {
      if (p %in% names(design_grid)) {
        env_max[[p]] <- max(design_grid[[p]], na.rm = TRUE)
      }
    }
    n_stage2_max_grid <- if (!is.na(n_s2_name)) {
      eval_resource(resources[[n_s2_name]], env_max)
    } else {
      max(n_stage2, na.rm = TRUE)
    }

    # Maximum stage 2 actually used by decision rules
    n_stage2_max_used <- max(n_stage2, na.rm = TRUE)

    # Cluster calculations
    clusters_s1 <- if (!is.na(clusters_s1_name)) {
      eval_resource(resources[[clusters_s1_name]], env_base)
    } else {
      n_arms * (env_base$k1 %||% 0)
    }

    # Expected stage 2 clusters
    clusters_s2_vec <- sapply(seq_len(nrow(opt_designs)), function(i) {
      if (!opt_designs$continue[i]) return(0)

      env <- env_base
      for (p in stage2_params) {
        if (p %in% names(opt_designs)) {
          val <- opt_designs[[p]][i]
          env[[p]] <- if (is.na(val)) 0 else val
        }
      }

      if (!is.na(clusters_s2_name)) {
        eval_resource(resources[[clusters_s2_name]], env)
      } else {
        n_arms * (env$k2 %||% 0)
      }
    })

    max_clusters_s2 <- if (!is.na(clusters_s2_name)) {
      eval_resource(resources[[clusters_s2_name]], env_max)
    } else {
      max(clusters_s2_vec, na.rm = TRUE)
    }

    # Cost calculations using generated functions
    cost_fn_s2 <- generate_cost_fn(spec, "stage2")
    cost_fn_total <- generate_cost_fn(spec, "total")

    # Stage 1 cost
    cost_stage1 <- cost_fn_total(model1, rho) - cost_fn_s2(model1, rho)

    # Stage 2 costs (vector)
    cost_stage2 <- sapply(seq_len(nrow(opt_designs)), function(i) {
      if (!opt_designs$continue[i]) return(0)

      env <- as.list(model1)
      for (p in stage2_params) {
        if (p %in% names(opt_designs)) {
          val <- opt_designs[[p]][i]
          env[[p]] <- if (is.na(val)) 0 else val
        }
      }

      cost_fn_s2(env, rho)
    })

    # Expected values using quadrature weights
    continue_idx <- opt_designs$continue
    weights <- results$quadrature$weights

    E_clusters_s2 <- sum(clusters_s2_vec[continue_idx] * weights[continue_idx])
    E_cost_stage2 <- sum(cost_stage2 * weights)
    E_cost_total <- cost_stage1 + E_cost_stage2

    # Maximum cost
    cost_stage2_max <- cost_fn_s2(env_max, rho)
    cost_total_max <- cost_stage1 + cost_stage2_max

    list(
      n_stage1 = n_stage1,
      n_stage2 = n_stage2,
      n_total_max = n_stage1 + n_stage2_max_used,      # Actual max from decisions
      n_total_max_grid = n_stage1 + n_stage2_max_grid, # Theoretical max from grid
      n_total_min = n_stage1,
      metrics = list(
        clusters_s1 = clusters_s1,
        max_clusters_s2 = max_clusters_s2,
        E_clusters_s2 = E_clusters_s2,
        E_clusters = clusters_s1 + E_clusters_s2,
        max_clusters = clusters_s1 + max_clusters_s2,
        cost_stage1 = cost_stage1,
        E_cost_stage2 = E_cost_stage2,
        E_cost = E_cost_total,
        max_cost = cost_stage1 + max(cost_stage2, na.rm = TRUE),
        max_cost_grid = cost_total_max  # Also add this
      )
    )
  }
}

# =============================================================================
# Model Factory Generation
# =============================================================================

#' Generate model function factory from specification
#'
#' @param spec A crt_design_spec object
#' @return A factory function that creates model functions with fresh caches
generate_model_factory <- function(spec) {

  model_builder <- spec$model_builder

  function() {
    cache <- new.env()
    cache$mod <- NULL
    cache$cache_key <- NULL

    function(design_params, fixed_params) {
      fixed_params$cache <- cache
      model_builder(design_params, fixed_params)
    }
  }
}

# =============================================================================
# Fixed Design Comparison
# =============================================================================

#' Generate fixed design model from adaptive spec
#'
#' @param spec A crt_design_spec object
#' @return A model builder function for single-stage designs
generate_fixed_model_builder <- function(spec) {

  # The fixed model builder combines stage 1 and stage 2 into one
  # This is design-specific, so we provide defaults for common types

  if (spec$design_type == "parallel") {
    return(fixed_parallel_model_builder)
  } else if (spec$design_type == "crossover") {
    return(fixed_crossover_model_builder)
  } else if (spec$design_type == "stepped_wedge") {
    return(fixed_sw_model_builder)
  }

  # Generic fallback - user must provide
  NULL
}

#' Generate fixed design cost function
#'
#' @param spec A crt_design_spec object
#' @return A cost function for single-stage designs
generate_fixed_cost_fn <- function(spec) {

  resources <- spec$resources
  cost_structure <- spec$cost_structure
  n_arms <- spec$n_arms

  function(model_output, rho) {
    env <- as.list(model_output)
    env$rho <- rho
    env$n_arms <- n_arms

    # For fixed designs, use n_total and clusters directly
    n_cost <- env$n_total %||% (n_arms * (env$k %||% env$k1 %||% 1) * (env$m %||% env$m1 %||% 1))
    k_cost <- n_arms * (env$k %||% env$k1 %||% 0)

    n_cost + rho * k_cost
  }
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
