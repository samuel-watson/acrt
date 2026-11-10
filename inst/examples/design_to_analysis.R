# =============================================================================
# From design to analysis: a simulated adaptive cluster randomised trial
# =============================================================================
# A two-arm parallel trial with a binary outcome, randomised within two regions.
# The script designs the trial, simulates stage 1, makes the interim decision,
# simulates the stage 2 that was chosen, and reports the final analysis.
#
# Run it with:
#   source(system.file("examples", "design_to_analysis.R", package = "acrt"))

library(acrt)
set.seed(2043)

# -----------------------------------------------------------------------------
# 1. Design the trial
# -----------------------------------------------------------------------------
# The stage 2 grid lists the options available at the interim: new clusters per
# arm (k2) and participants per cluster-period (m2).

design <- parallel_crt(
  icc = 0.05,
  cac = 0.8,
  delta = -0.10,          # planned risk difference
  baseline = 0.30,
  k1 = 16,                # stage 1 clusters per arm
  m1 = 30,                # stage 1 participants per cluster
  k2 = 0:4,
  m2 = c(20, 30, 40, 50),
  rho = 30,
  family = "binomial"
)

planned <- adaptive_analysis(design, target_power = 0.8, tol = 0.01,
                             method = "lambda", verbose = FALSE)
print(planned)

# -----------------------------------------------------------------------------
# 2. Fix the analysis plan
# -----------------------------------------------------------------------------
# The family, weights, boundaries and planned autocorrelation come from the
# design. For a parallel_crt() design the schedule has four sequences:
#   1 = control, 2 = intervention (stage 1 clusters)
#   3 = control, 4 = intervention (clusters added at stage 2)

plan <- analysis_plan(planned, strata = "region")
print(plan)

# -----------------------------------------------------------------------------
# 3. Randomise and simulate stage 1
# -----------------------------------------------------------------------------
# The simulation uses the design's own model: a cluster random effect plus a
# cluster-period random effect on the logit scale.

truth <- convert_to_glm_scale("binomial", baseline = 0.30, delta = -0.10,
                              icc = 0.05, cac = 0.8)
tau_cluster <- sqrt(truth$cov_pars[1])
tau_period <- sqrt(truth$cov_pars[2])
region_effect <- c(North = -0.15, South = 0.15)

# Randomise 16 clusters to each arm, balanced within region
alloc <- data.frame(
  cluster = 1:32,
  region = rep(c("North", "South"), each = 16),
  sequence = c(sample(rep(1:2, 8)), sample(rep(1:2, 8)))
)
cluster_effect <- rnorm(32, 0, tau_cluster)

simulate_period <- function(alloc, cluster_effect, period, m) {
  on_intervention <- plan$schedule[alloc$sequence, period]
  logit <- truth$beta[1] + region_effect[alloc$region] +
    truth$beta[2] * on_intervention + cluster_effect[alloc$cluster] +
    rnorm(nrow(alloc), 0, tau_period)
  out <- data.frame(
    cluster = rep(alloc$cluster, each = m),
    period = period,
    outcome = rbinom(nrow(alloc) * m, 1, rep(plogis(logit), each = m))
  )
  # 5% of participants have no recorded primary outcome
  out$outcome[runif(nrow(out)) < 0.05] <- NA
  out
}

data_stage1 <- simulate_period(alloc, cluster_effect, period = 1,
                               m = design$stage1_grid$m1)

# -----------------------------------------------------------------------------
# 4. Interim analysis
# -----------------------------------------------------------------------------
# interim_statistic() locks the nuisance parameters. interim_decision() passes
# z1 and the updated ICC, autocorrelation and baseline to interim_analysis().

interim <- interim_statistic(plan, data_stage1, alloc)
print(interim)

decision <- interim_decision(interim, verbose = FALSE)
print(decision)

# -----------------------------------------------------------------------------
# 5. Stage 2 and final analysis, or stop at the interim
# -----------------------------------------------------------------------------

if (decision$decision == "continue") {

  k2 <- decision$stage2_design$k2
  m2 <- decision$stage2_design$m2

  # Clusters added at stage 2, k2 per arm, spread across regions
  if (k2 > 0) {
    new_clusters <- data.frame(
      cluster = 32 + seq_len(2 * k2),
      region = rep_len(c("North", "South"), 2 * k2),
      sequence = rep(3:4, each = k2)
    )
    alloc_final <- rbind(alloc, new_clusters)
    cluster_effect <- c(cluster_effect, rnorm(2 * k2, 0, tau_cluster))
  } else {
    alloc_final <- alloc
  }

  data_stage2 <- simulate_period(alloc_final, cluster_effect, period = 2, m = m2)
  data_final <- rbind(data_stage1, data_stage2)

  final <- final_statistic(interim, data_final, alloc = alloc_final)
  print(final)

  estimates <- effect_estimates(final)

} else {

  # With an efficacy or futility stop, the analysis uses stage 1 data only
  estimates <- effect_estimates(interim)
}

# -----------------------------------------------------------------------------
# 6. Results
# -----------------------------------------------------------------------------
# The `adaptive` row is the primary result. Its interval is the repeated
# confidence interval. The other rows are descriptive.

print(estimates)
cat(sprintf("\nTrue log odds ratio used in the simulation: %.3f\n", truth$beta[2]))

# Stage-wise statistics at any hypothesised effect, on the log odds ratio scale
if (exists("final")) print(unlist(stage_statistics(final, beta = -0.25)))

# -----------------------------------------------------------------------------
# Other designs
# -----------------------------------------------------------------------------
# For a design built with crt_design_spec(), give the schedule and the stage 1
# periods. Use NA for periods in which a sequence is not part of the trial.
# A stepped-wedge design has period effects added automatically, because
# treatment is confounded with calendar time:
#
#   schedule <- rbind(c(0, 1, 1, 1, 1),
#                     c(0, 0, 1, 1, 1),
#                     c(0, 0, 0, 1, 1),
#                     c(0, 0, 0, 0, 1))
#   plan <- analysis_plan(planned, schedule = schedule, stage1 = 1:3)
#
# If the interim decision changes the stage 2 rollout or the number of periods,
# pass the realised schedule to final_statistic(interim, data, schedule = ...).
