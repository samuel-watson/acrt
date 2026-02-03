# =============================================================================
# Adaptive Cluster Trials Package
# =============================================================================
# A framework for designing and analysing adaptive cluster randomised trials
# using score test decomposition.
#
# Main user-facing functions:
#   - parallel_crt(), crossover_crt(), stepped_wedge_crt(): Create design objects
#   - adaptive_analysis(): Run the analysis
#   - find_pareto(): Find Pareto-optimal designs
#   - compare_to_fixed(): Compare to fixed design
#
# Quick start:
#   design <- parallel_crt(icc = 0.05, delta = 0.25, k1 = 8:12, m1 = seq(10, 50, 10))
#   results <- adaptive_analysis(design, target_power = 0.8)
#   plot(results)

# =============================================================================
# Load Components
# =============================================================================

# Source order matters - dependencies first

# 1. Core computational functions (no dependencies)
# source("R/adaptive_cluster_trial_v2.R")
#
# # 2. Design specification system
# source("R/design_spec.R")
#
# # 3. Design templates (depends on design_spec.R)
# source("R/design_templates.R")
#
# # 4. Design space exploration and plotting
# source("R/adaptive_design_space_v2.R")
#
# # 5. Simplified analysis interface (depends on all above)
# source("R/analysis_interface.R")

# =============================================================================
# Package Information
# =============================================================================

.pkg_info <- list(
  name = "adaptiveCRT",
  version = "0.2.0",
  description = "Adaptive Cluster Randomised Trial Design",
  author = "Sam Watson",
  depends = c("glmmrBase", "Matrix", "statmod", "ggplot2", "scales")
)

#' Print package information
#' @export
adaptive_crt_info <- function() {
  cat(sprintf("%s v%s\n", .pkg_info$name, .pkg_info$version))
  cat(.pkg_info$description, "\n\n")
  cat("Available design constructors:\n")
  cat("  parallel_crt()      - Parallel cluster trial\n")
  cat("  crossover_crt()     - Crossover cluster trial\n")
  cat("  stepped_wedge_crt() - Stepped-wedge (rollout) trial\n\n")
  cat("Main functions:\n")
  cat("  adaptive_analysis() - Run analysis\n")
  cat("  find_pareto()       - Find Pareto frontier\n")
  cat("  compare_to_fixed()  - Compare to fixed design\n")
}

# =============================================================================
# Convenience Re-exports
# =============================================================================

# These are the main user-facing functions, re-exported for discoverability

#' @rdname parallel_crt
#' @export
parallel <- parallel_crt

#' @rdname crossover_crt
#' @export
crossover <- crossover_crt

#' @rdname stepped_wedge_crt
#' @export
stepped_wedge <- stepped_wedge_crt

#' @rdname stepped_wedge_crt
#' @export
rollout <- stepped_wedge_crt

# =============================================================================
# Startup Message
# =============================================================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("adaptiveCRT v%s loaded", .pkg_info$version))
  packageStartupMessage("Use adaptive_crt_info() for help")
}

cat("adaptiveCRT framework loaded\n")
cat("Use adaptive_crt_info() for available functions\n")
