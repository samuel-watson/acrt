# =============================================================================
# Example: Simplified User Interface
# =============================================================================
# This example demonstrates how easy it is to use the framework with the
# new simplified interface. Compare this to the 150+ lines previously required!

# Load the framework
setwd("~/adaptive_trials")  # Adjust path as needed
source("adaptiveCRT.R")

# =============================================================================
# Example 1: Basic Parallel CRT (Simplest Usage)
# =============================================================================

cat("\n=== Example 1: Basic Parallel CRT ===\n\n")

# Define the design in ONE call - no custom functions needed!
design <- parallel_crt(
  icc = 0.05,           # Intra-cluster correlation
  cac = 0.8,            # Cluster autocorrelation
  delta = 0.25,         # Treatment effect
  k1 = c(8, 10, 12, 14, 16),    # Stage 1 clusters per arm (explore these)
  m1 = seq(10, 50, 5), # Stage 1 individuals per cluster
  k2 = 0:4,             # Stage 2 new clusters (search space)
  rho = 30              # Cluster-to-individual cost ratio
)

# Run the analysis - ONE function call
results <- adaptive_analysis(design, target_power = 0.8, verbose = TRUE)

# View summary
summary(results)

# Plot results
plot(results, type = "EN")      # Expected sample size
plot(results, type = "Nmax")    # Maximum sample size
plot(results, type = "pareto")  # Pareto frontier

# =============================================================================
# Example 2: Single Design Analysis
# =============================================================================

cat("\n=== Example 2: Single Design Analysis ===\n\n")

# For a single specific design, just provide scalars
design_single <- parallel_crt(
  icc = 0.05,
  delta = 0.25,
  k1 = 14,    # Single value
  m1 = 15,    # Single value
  k2 = 0:6,
  rho = 30
)

# Analyse
results_single <- adaptive_analysis(design_single, target_power = 0.8)

# Summary shows more detail for single designs
summary(results_single)

# Decision rules
plot(results_single, type = "decision")

# Get decision rules as data frame
rules <- get_decision_rules(results_single)
print(rules[seq(1, nrow(rules), by = 5), c("z1", "continue", "cp", "m2", "k2")])

# =============================================================================
# Example 3: Pareto Frontier Analysis
# =============================================================================

cat("\n=== Example 3: Pareto Frontier ===\n\n")

# Find Pareto-optimal designs (minimise E[N] and N_max)
pareto <- find_pareto(results, objectives = list(E_N = "min", N_max = "min"))
print(pareto[, c("k1", "m1", "E_N", "N_max", "power_stage1")])

# Different objectives: minimise cost, maximise early stopping
pareto2 <- find_pareto(results,
                       objectives = list(E_cost = "min", prob_efficacy = "max"))
print(pareto2[, c("k1", "m1", "E_cost", "prob_efficacy")])

# With constraints
pareto3 <- find_pareto(results,
                       objectives = list(E_N = "min", N_max = "min"),
                       constraints = list(power_stage1 = c(0.3, 0.6)))
print(pareto3[, c("k1", "m1", "E_N", "N_max", "power_stage1")])

# =============================================================================
# Example 4: Compare to Fixed Design
# =============================================================================

cat("\n=== Example 4: Compare to Fixed Design ===\n\n")

# Compare adaptive designs to optimal fixed design

comparison <- compare_to_fixed(results)
print(comparison)

# Plot comparison
plot(comparison, type = "savings")
plot(comparison, type = "frontier")

# =============================================================================
# Example 5: Stepped-Wedge Design
# =============================================================================

cat("\n=== Example 5: Stepped-Wedge Design ===\n\n")

# Stepped-wedge with rollout adaptation
sw_design <- stepped_wedge_crt(
  icc = 0.05,
  cac = 0.8,
  delta = 0.3,
  K = 24,           # Total clusters (fixed)
  t1 = 3,           # Stage 1 periods
  r1 = 4,           # Clusters switching per period in stage 1
  m = 20,           # Individuals per cluster-period
  t2 = 1:5,         # Additional periods to search over
  r2 = c(2, 4, 6),  # Rollout rates to search over
  rho = 100         # Period-to-individual cost ratio
)

# Analyse
sw_results <- adaptive_analysis(sw_design, target_power = 0.8, verbose = TRUE)
summary(sw_results)

# =============================================================================
# Example 6: Crossover Design
# =============================================================================

cat("\n=== Example 6: Crossover Design ===\n\n")

xo_design <- crossover_crt(
  icc = 0.05,
  cac = 0.8,
  delta = 0.25,
  k1 = c(6, 8, 10),
  m1 = seq(20, 60, 10),
  k2 = 0:3,
  rho = 30
)

xo_results <- adaptive_analysis(xo_design, target_power = 0.8, verbose = TRUE)
summary(xo_results)
plot(xo_results, type = "EN")

# =============================================================================
# Example 7: Custom Design (Advanced)
# =============================================================================

cat("\n=== Example 7: Custom Design (Advanced) ===\n\n")

# For non-standard designs, create a custom specification
# This still avoids most boilerplate

custom_spec <- crt_design_spec(
  stage1_params = c("k1", "m1"),
  stage2_params = c("k2", "m2"),
  resources = list(
    n_s1 = ~ 2 * k1 * m1,
    n_s2 = ~ 2 * k1 * m2 + 2 * k2 * m2,
    clusters_s1 = ~ 2 * k1,
    clusters_s2 = ~ 2 * k2
  ),
  cost_structure = list(
    weights = c(n = 1, clusters = "rho"),
    stage2_resources = c("n_s2", "clusters_s2")
  ),
  model_builder = parallel_model_builder,  # Can use built-in or custom
  design_type = "parallel"
)

# Create design manually
custom_design <- structure(
  list(
    spec = custom_spec,
    fixed_params = list(icc = 0.05, cac = 0.8, delta = 0.25),
    stage1_grid = expand.grid(k1 = 8:12, m1 = c(20, 30, 40)),
    stage2_grid_fn = function(s1) expand.grid(m2 = seq(10, 60, 10), k2 = 0:4),
    rho = 30
  ),
  class = c("custom_crt", "crt_design")
)

custom_results <- adaptive_analysis(custom_design, target_power = 0.8)
summary(custom_results)

# =============================================================================
# Summary: Lines of Code Comparison
# =============================================================================

cat("\n=== Code Comparison ===\n\n")
cat("OLD approach (per design type):\n")
cat("  - parallel_model(): ~80 lines\n")
cat("  - make_parallel_model_fn(): ~10 lines\n")
cat("  - parallel_cost_fn(): ~10 lines\n")
cat("  - parallel_sample_size_fn(): ~50 lines\n")
cat("  - fixed_parallel_model(): ~40 lines\n")
cat("  - fixed_parallel_cost_fn(): ~5 lines\n")
cat("  - Analysis code: ~30 lines\n")
cat("  TOTAL: ~225 lines per design type\n\n")

cat("NEW approach:\n")
cat("  design <- parallel_crt(icc=0.05, delta=0.25, k1=8:12, m1=10:50, rho=30)\n")
cat("  results <- adaptive_analysis(design, target_power = 0.8)\n")
cat("  TOTAL: 2 lines!\n\n")

cat("Reduction: 99%\n")
