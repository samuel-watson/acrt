# =============================================================================
# Package documentation, imports and start-up
# =============================================================================

#' acrt: Adaptive Cluster Randomised Trial Design and Analysis
#'
#' Tools for designing and analysing two-stage adaptive cluster randomised
#' trials using a combination test of stage-wise efficient scores.
#'
#' Designs are created with [parallel_crt()] or, for any other layout, with
#' [crt_design_spec()]. [adaptive_analysis()] calibrates the decision rule and
#' [interim_analysis()] applies it during a trial. [analysis_plan()],
#' [interim_statistic()], [final_statistic()] and [effect_estimates()] analyse
#' the trial data.
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom stats aggregate as.formula binomial complete.cases dnorm dt
#'   gaussian integrate optimize pnorm poisson pt qnorm qt uniroot
#' @importFrom utils globalVariables
"_PACKAGE"

# Column names used inside ggplot2 aes() and subset() calls.
utils::globalVariables(c(
  "E_cost", "E_cost_ratio", "E_savings_pct", "cost", "cp", "cum_prob",
  "density_scaled", "design_scaled", "hypothesis", "line_type", "max_cost",
  "max_cost_ratio", "max_within_fixed", "metric", "n_total", "outcome",
  "p_reject", "power_stage1", "prob_continue", "prob_efficacy", "prob_futility",
  "probability", "scenario", "value", "weight", "ymax", "ymin", "z1"
))

#' @rdname parallel_crt
#' @export
parallel <- parallel_crt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf("%s %s: see vignette(package = \"%s\")", pkgname,
            utils::packageVersion(pkgname), pkgname)
  )
}
