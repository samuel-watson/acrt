# Adaptive cluster randomised trials
This R package provide functionality to analyse two-stage adaptive cluster randomised trial designs. Full details of the methodology is provided in Watson & Martin (2026). The package allows for specification of any set of trial designs with arbitrary correlation structure. We provide in-built functionality for parallel trials and incorporate two complete vignettes demonstrating how to use the more flexible aspects of the package. The package is not currently submitted to CRAN while further tests and functions are added. To install please use
```
devtools::install_github("samuel-watson/acrt")
```
The vignettes can be loaded as:
```
vignette('adaptive_stepped_wedge', package = 'acrt')
vignette('adaptive_crossover', package = 'acrt')
```

## Parallel trial example
The code below reproduces the parallel trial example from the paper:
```
design <- parallel_crt(
  icc = 0.05,           # Intra-cluster correlation
  cac = 0.8,            # Cluster autocorrelation
  delta = 0.25,         # Treatment effect
  k1 = c(12:16),        # Stage 1 clusters per arm
  m1 = seq(10, 50, 5),  # Stage 1 individuals per cluster
  k2 = 0:4,             # Stage 2 new clusters (search space)
  rho = 30              # Cluster-to-individual cost ratio
)

# max cost minimisation
results <- adaptive_analysis(design, target_power = 0.8, verbose = TRUE, tol = 0.01, method = "cost_cap")
summary(results)

# Plot results
p1 <- plot(results, type = "EN")      # Expected sample size
p2 <- plot(results, type = "Nmax")    # Maximum sample size
p3 <- plot(results, type = "pareto", objectives = list(E_cost = "min", max_cost = "min") )  # Pareto frontier

pareto1 <- find_pareto(results, objectives = list(E_cost = "min", max_cost = "min"))
print(pareto1[, c("k1", "m1", "E_N", "N_max", "power_stage1")])

# expected cost minimisation
results2 <- adaptive_analysis(design, target_power = 0.8, verbose = TRUE, tol = 0.01, method = "lambda")
summary(results2)

# Plot results
p4 <- plot(results2, type = "EN")      # Expected sample size
p5 <- plot(results2, type = "Nmax")    # Maximum sample size
p6 <- plot(results2, type = "pareto", objectives = list(E_cost = "min", max_cost = "min") )  # Pareto frontier

pareto2 <- find_pareto(results, objectives = list(E_cost = "min", max_cost = "min"))
print(pareto2[, c("k1", "m1", "E_N", "N_max", "power_stage1")])

require(patchwork)

layout <- "
abc
def
"

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = "collect", design = layout) + plot_annotation(theme = theme(legend.position = "bottom")) & labs(title = NULL, subtitle = NULL)

# examine the chosen design2
design_single <- parallel_crt(
  icc = 0.05,
  delta = 0.25,
  k1 = 15,
  m1 = 20,
  k2 = 0:4,
  rho = 30
)

design_single2 <- parallel_crt(
  icc = 0.05,
  delta = 0.25,
  k1 = 16,
  m1 = 25,
  k2 = 0:4,
  rho = 30
)

# Analyse
results_single <- adaptive_analysis(design_single, target_power = 0.8, tol = 0.01, method = "lambda")
results_single2 <- adaptive_analysis(design_single2, target_power = 0.8, tol = 0.01, method = "cost_cap")

# Decision rules
q1 <- plot(results_single, type = "decision")
q2 <- plot(results_single2, type = "decision")

q1 + q2

# Get decision rules as data frame
rules <- get_decision_rules(results_single)
print(rules[seq(1, nrow(rules), by = 5), c("z1", "continue", "cp", "m2", "k2")])
```
