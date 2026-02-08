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
  k1 = c(12:16),    # Stage 1 clusters per arm (explore these)
  m1 = seq(10, 50, 5), # Stage 1 individuals per cluster
  k2 = 0:4,             # Stage 2 new clusters (search space)
  rho = 30              # Cluster-to-individual cost ratio
)

# Run the analysis - ONE function call
results <- adaptive_analysis(design, target_power = 0.8,
                             verbose = TRUE, tol = 0.01, method = "cost_cap")

# View summary
summary(results)

# Plot results
p1 <- plot(results, type = "EN")      # Expected sample size
p2 <- plot(results, type = "Nmax")    # Maximum sample size
p3 <- plot(results, type = "pareto", objectives = list(E_cost = "min", max_cost = "min") )  # Pareto frontier


pareto <- find_pareto(results, objectives = list(E_cost = "min", max_cost = "min"))

# For a single specific design, just provide scalars
design_single <- parallel_crt(
  icc = 0.05,
  delta = 0.25,
  k1 = 15,    # Single value
  m1 = 20,    # Single value
  k2 = 0:4,
  rho = 30
)

# Analyse
results_single2 <- adaptive_analysis(design_single, target_power = 0.8, tol = 0.01)

plot_cost_distribution(results_single2,results_cap = results_single,fixed_cost = comparison$fixed$optimal$cost,
                       sample_size_fn = results_single$sample_size_fn,rho = 30)

# Summary shows more detail for single designs
summary(results_single)

# Decision rules
p4b <- plot(results_single, type = "decision")

# Get decision rules as data frame
rules <- get_decision_rules(results_single)
print(rules[seq(1, nrow(rules), by = 5), c("z1", "continue", "cp", "m2", "k2")])



comparison <- compare_to_fixed(results)
print(comparison)

# Plot comparison
plot(comparison, type = "savings")
plot(comparison, type = "frontier")
p3 <- plot(comparison, type = "cost_ratio")


require(patchwork)

layo <- "
abc
def
"

p1 + p2 + p3 + p1b + p2b + p3b + plot_layout(design = layo, guides = "collect") &
  labs(title = NULL, subtitle = NULL) & theme(legend.position = "bottom")

p4 + p4b + plot_layout(ncol = 2, nrow = 6)
