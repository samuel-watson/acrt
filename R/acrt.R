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
  name = "acrt",
  version = "0.2.0",
  description = "Adaptive Cluster Randomised Trial Design",
  author = "Sam Watson",
  depends = c("glmmrBase", "Matrix", "statmod", "ggplot2", "scales")
)

# =============================================================================
# Convenience Re-exports
# =============================================================================

# These are the main user-facing functions, re-exported for discoverability

#' @rdname parallel_crt
#' @export
parallel <- parallel_crt

# =============================================================================
# Startup Message
# =============================================================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("adaptiveCRT v%s loaded", .pkg_info$version))
  packageStartupMessage("Vignettes demonstrate use of the package:\nvignette('adaptive_stepped_wedge', package = 'acrt')\nvignette('adaptive_crossover', package = 'acrt')")
}

cat("adaptive CRT framework loaded\n")
