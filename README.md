# Custom Design Guide for adaptiveCRT

This guide explains how to create custom adaptive trial designs using the adaptiveCRT framework. You only need custom code if the built-in designs (`parallel_crt`, `crossover_crt`, `stepped_wedge_crt`) don't fit your needs.

**For quick reference in R, run:** `custom_design_help()`

## Overview

To create a custom design, you need to provide:

| Component | Required? | Description |
|-----------|-----------|-------------|
| `model_builder` | **Yes** | Computes efficient information decomposition |
| `resources` | **Yes** | Formulas for sample sizes, clusters, etc. |
| `cost_fn` | No | Auto-generated from resources by default |
| `sample_size_fn` | No | Auto-generated from resources by default |

## Quick Start Template

```r
# 1. Define your model builder (see detailed docs below)
my_model_builder <- function(design_params, fixed_params) {
  # ... your code here ...
  list(I1_eff=..., I2_eff=..., I_eff=..., w1=..., w2=..., b1=...)
}

# 2. Create design specification
my_spec <- crt_design_spec(
  stage1_params = c("k1", "m1"),
  stage2_params = c("k2", "m2"),
  resources = list(
    n_s1 = ~ 2 * k1 * m1,
    n_s2 = ~ 2 * k1 * m2 + 2 * k2 * m2,
    clusters_s1 = ~ 2 * k1,
    clusters_s2 = ~ 2 * k2
  ),
  model_builder = my_model_builder,
  design_type = "my_custom_type"
)

# 3. Create design object
my_design <- structure(
  list(
    spec = my_spec,
    fixed_params = list(icc = 0.05, cac = 0.8, delta = 0.25),
    stage1_grid = expand.grid(k1 = 8:12, m1 = c(20, 30, 40)),
    stage2_grid_fn = function(s1) expand.grid(m2 = 10:50, k2 = 0:4),
    rho = 30
  ),
  class = c("my_custom_crt", "crt_design")
)

# 4. Analyse
results <- adaptive_analysis(my_design, target_power = 0.8)
```

---

## Model Builder Function

The `model_builder` is the core function that computes the statistical quantities needed for adaptive decision-making.

### Signature

```r
model_builder <- function(design_params, fixed_params) { ... }
```

### Inputs

#### `design_params`
A single-row `data.frame` with stage 2 parameters from your stage 2 grid.

```r
# Example
design_params <- data.frame(m2 = 30, k2 = 2)
design_params$m2  # Access values like this
```

#### `fixed_params`
A list containing:

| Field | Source | Description |
|-------|--------|-------------|
| Stage 1 params | `stage1_grid` | e.g., `k1`, `m1`, `K`, `t1` |
| Constants | `fixed_params_base` | e.g., `icc`, `cac`, `delta` |
| `cache` | Framework | Environment for model caching |

```r
fixed_params$k1      # Stage 1 clusters
fixed_params$icc     # Intra-cluster correlation
fixed_params$delta   # Treatment effect
fixed_params$cache   # Model cache environment
```

### Required Outputs

The function **must** return a list with these fields:

| Field | Type | Description |
|-------|------|-------------|
| `I1_eff` | scalar | Stage 1 efficient information |
| `I2_eff` | scalar | Stage 2 efficient information |
| `I_eff` | scalar | Total efficient information |
| `w1` | scalar | Stage 1 weight = `sqrt(I1_eff / I_eff)` |
| `w2` | scalar | Stage 2 weight = `sqrt(I2_eff / I_eff)` |
| `b1` | scalar | Treatment effect (same as `delta`) |

**Mathematical interpretation:**
- Stage 1 score: Z₁ ~ N(b₁ × √I₁_eff, 1) under H₁
- Combination test: Z = w₁Z₁ + w₂Z₂|₁
- Weights satisfy: w₁² + w₂² = 1

### Optional Outputs

Include any quantities needed by cost or sample size calculations:

```r
list(
  # ... required fields ...
  
  # Sample size tracking
  n1 = 480,        # Observations in stage 1
  n2 = 320,        # Observations in stage 2
  
  # Design parameters (for resource formulas)
  k1 = 8, k2 = 2,
  m1 = 30, m2 = 20,
  
  # Any other needed quantities
  T_total = 6
)
```

### Computing Efficient Information

Use the `efficient_score_decomposition()` helper:

```r
# Build your design matrices
X <- cache$mod$mean$X          # Fixed effects design matrix
V <- compute_V_matrix(...)     # Marginal covariance matrix

# Define stage indices
idx1 <- 1:n1                   # Row indices for stage 1
idx2 <- (n1 + 1):(n1 + n2)     # Row indices for stage 2

# j = column index of treatment effect in X
eff <- efficient_score_decomposition(X, V, idx1, idx2, j = 2)

# Use results
I1_eff <- eff$I1_eff
I2_eff <- eff$I2_eff
I_eff <- eff$I_eff
```

The "efficient" information projects out nuisance parameters (intercepts, period effects) that aren't fully estimable within each stage alone.

### Model Caching (Critical for Performance)

The `model_builder` is called many times. Building glmmrBase models is slow. **Always use caching:**

```r
cache <- fixed_params$cache
cache_key <- paste(k1, k2, sep = "_")  # Unique ID for model structure

if (is.null(cache$mod) || cache$cache_key != cache_key) {
  # Model structure changed - must rebuild
  cache$mod <- Model$new(
    formula = ~ trt + (1|gr(cl)),
    data = df,
    family = gaussian(),
    mean = c(0, delta),
    covariance = icc,
    weights = df$n,
    var_par = 1 - icc
  )
  cache$cache_key <- cache_key
} else {
  # Only weights/parameters changed - just update (fast)
  glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
  cache$mod$update_parameters(cov.pars = new_variance_params)
}
```

**Rule:** Rebuild when cluster structure changes; update when only sample sizes change.

### Computing V Matrix from glmmrBase

Standard pattern:

```r
X <- cache$mod$mean$X           # Fixed effects design matrix
D <- cache$mod$covariance$D     # Random effects covariance
Z <- cache$mod$covariance$Z     # Random effects design matrix

library(Matrix)
Z_sp <- Matrix(Z, sparse = TRUE)
D_sp <- Matrix(D, sparse = TRUE)

# Marginal covariance: V = σ²I + ZDZ'
S <- Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% tcrossprod(D_sp, Z_sp)
V <- as.matrix(S)
```

---

## Resources Specification

Resources are formulas that define how to calculate quantities like sample sizes, cluster counts, and time periods.

### Naming Convention

```
{resource_type}_{stage}
```

| Name | Description |
|------|-------------|
| `n_s1` | Observations in stage 1 |
| `n_s2` | Observations in stage 2 |
| `clusters_s1` | Clusters in stage 1 |
| `clusters_s2` | **New** clusters in stage 2 |
| `periods_s2` | Additional time periods in stage 2 |

### Formula Syntax

Use R formula notation:

```r
resources = list(
  n_s1        = ~ n_arms * k1 * m1,
  n_s2        = ~ n_arms * k1 * m2 + n_arms * k2 * m2,
  clusters_s1 = ~ n_arms * k1,
  clusters_s2 = ~ n_arms * k2
)
```

### Available Variables

- All stage 1 parameters (`k1`, `m1`, `K`, `t1`, etc.)
- All stage 2 parameters (`k2`, `m2`, `t2`, `r2`, etc.)
- `n_arms` (from design spec)
- `rho` (cost ratio, when evaluating costs)
- Any outputs from `model_builder`

### Cost Structure

```r
cost_structure = list(
  weights = c(n = 1, clusters = "rho", periods = "rho"),
  stage2_resources = c("n_s2", "clusters_s2")
)
```

| Field | Description |
|-------|-------------|
| `weights` | Maps resource types to cost multipliers |
| `stage2_resources` | Which resources to include for λ-calibration |

**Weight values:**
- Numeric (e.g., `1`): used directly
- Character (e.g., `"rho"`): looked up as variable

---

## Custom Cost Function (Optional)

Only needed if the auto-generated cost doesn't fit your design.

### Signature

```r
cost_fn <- function(model_output, rho) { ... }
```

### Inputs

| Parameter | Description |
|-----------|-------------|
| `model_output` | List returned by your `model_builder` |
| `rho` | Cost ratio from `design$rho` |

### Output

Return a **scalar**: the **stage 2 cost only** (not total cost).

This is used in the criterion: CP(z₁) - λ × cost

### Example

```r
my_cost_fn <- function(model_output, rho) {
  # Observations cost 1 each
  n_cost <- model_output$n2
  
  # Clusters cost rho each (2 arms)
  k_cost <- 2 * model_output$k2
  
  # Return stage 2 cost
  n_cost + rho * k_cost
}
```

---

## Custom Sample Size Function (Optional)

Only needed if the auto-generated sample size doesn't fit your design.

### Signature

```r
sample_size_fn <- function(opt_designs, results, rho = 30) { ... }
```

### Inputs

#### `opt_designs`
Data frame with one row per quadrature point (z₁ value):

| Column | Description |
|--------|-------------|
| `z1` | Stage 1 test statistic value |
| `continue` | TRUE if trial continues |
| `efficacy_stop` | TRUE if stopping for efficacy |
| `futility_stop` | TRUE if stopping for futility |
| `cp` | Conditional power |
| Stage 2 params | `m2`, `k2`, etc. (NA when not continuing) |

#### `results`
Full results object from `find_lambda_for_power`:

```r
results$models$list[[1]]       # First model (has stage 1 params)
results$models$design_grid     # Stage 2 design grid
results$quadrature$weights     # For computing expectations
```

### Required Output

```r
list(
  n_stage1    = scalar,                          # Stage 1 observations
  n_stage2    = vector,                          # Per z1 point
  n_total_max = scalar,                          # Maximum possible
  n_total_min = scalar                           # Minimum (usually = n_stage1)
)
```

### Optional Output

```r
metrics = list(
  E_cost = ...,
  max_cost = ...,
  E_clusters = ...
)
```

These appear in the summary table.

### Example

```r
my_sample_size_fn <- function(opt_designs, results, rho = 30) {
  
  model1 <- results$models$list[[1]]
  k1 <- model1$k1
  m1 <- model1$m1
  n_arms <- 2
  
  # Stage 1 (fixed)
  n_stage1 <- n_arms * k1 * m1
  
  # Stage 2 (varies by decision)
  m2 <- ifelse(is.na(opt_designs$m2), 0, opt_designs$m2)
  k2 <- ifelse(is.na(opt_designs$k2), 0, opt_designs$k2)
  
  n_stage2 <- ifelse(
    opt_designs$continue,
    n_arms * k1 * m2 + n_arms * k2 * m2,
    0
  )
  
  # Maximum from grid
  grid <- results$models$design_grid
  n_stage2_max <- n_arms * k1 * max(grid$m2) + n_arms * max(grid$k2) * max(grid$m2)
  
  # Costs
  weights <- results$quadrature$weights
  cost_stage1 <- n_stage1 + rho * (n_arms * k1)
  cost_stage2 <- ifelse(
    opt_designs$continue,
    n_stage2 + rho * (n_arms * k2),
    0
  )
  E_cost <- cost_stage1 + sum(cost_stage2 * weights)
  
  list(
    n_stage1 = n_stage1,
    n_stage2 = n_stage2,
    n_total_max = n_stage1 + n_stage2_max,
    n_total_min = n_stage1,
    metrics = list(
      E_cost = E_cost,
      max_cost = cost_stage1 + n_stage2_max + rho * (n_arms * max(grid$k2))
    )
  )
}
```

---

## Complete Example: Staggered Enrollment Design

A design where enrollment is staggered over time, with the option to extend enrollment in stage 2.

### Design Description

- K clusters total (fixed)
- Stage 1: Enroll clusters over t₁ periods at rate r₁ per period
- Stage 2: Optionally extend by t₂ periods at rate r₂
- Outcome measured at fixed time after enrollment

### Parameters

| Type | Parameters |
|------|------------|
| Stage 1 (fixed) | `K`, `t1`, `r1`, `m` |
| Stage 2 (adaptive) | `t2`, `r2` |

### Implementation

```r
# Step 1: Model builder
staggered_model_builder <- function(design_params, fixed_params) {
  
  K <- fixed_params$K
  t1 <- fixed_params$t1
  r1 <- fixed_params$r1
  m <- fixed_params$m
  t2 <- design_params$t2
  r2 <- design_params$r2
  icc <- fixed_params$icc
  delta <- fixed_params$delta
  cache <- fixed_params$cache
  
  T_total <- t1 + t2
  
  # Calculate enrollment times
  enroll_time <- numeric(K)
  enrolled <- 0
  
  for (t in 1:t1) {
    n_new <- min(r1, K - enrolled)
    if (n_new > 0) {
      enroll_time[(enrolled + 1):(enrolled + n_new)] <- t
      enrolled <- enrolled + n_new
    }
  }
  
  if (t2 > 0 && enrolled < K) {
    for (t in (t1 + 1):T_total) {
      n_new <- min(r2, K - enrolled)
      if (n_new > 0) {
        enroll_time[(enrolled + 1):(enrolled + n_new)] <- t
        enrolled <- enrolled + n_new
      }
      if (enrolled >= K) break
    }
  }
  
  enroll_time[enroll_time == 0] <- T_total + 1
  
  # Build data
  df <- data.frame(
    cl = 1:K,
    t = enroll_time,
    trt = rep(0:1, length.out = K)
  )
  df <- df[df$t <= T_total, ]
  df$n <- m
  
  idx1 <- which(df$t <= t1)
  idx2 <- which(df$t > t1)
  
  # Build model with caching
  cache_key <- paste(K, T_total, r2, sep = "_")
  
  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    cache$mod <- Model$new(
      ~ trt + (1|gr(cl)),
      data = df,
      family = gaussian(),
      mean = c(0, delta),
      covariance = icc,
      weights = df$n,
      var_par = 1 - icc
    )
    cache$cache_key <- cache_key
  } else {
    glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
  }
  
  # Compute V and efficient information
  X <- cache$mod$mean$X
  Z <- cache$mod$covariance$Z
  D <- cache$mod$covariance$D
  
  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) + 
       Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)
  
  eff <- efficient_score_decomposition(X, V, idx1, idx2, j = 2)
  
  list(
    I1_eff = eff$I1_eff,
    I2_eff = eff$I2_eff,
    I_eff = eff$I_eff,
    w1 = sqrt(eff$I1_eff / eff$I_eff),
    w2 = sqrt(eff$I2_eff / eff$I_eff),
    b1 = delta,
    n1 = length(idx1) * m,
    n2 = length(idx2) * m,
    clusters_s1 = length(idx1),
    clusters_s2 = length(idx2),
    K = K, t1 = t1, t2 = t2, m = m, T_total = T_total
  )
}

# Step 2: Create specification
staggered_spec <- crt_design_spec(
  stage1_params = c("K", "t1", "r1", "m"),
  stage2_params = c("t2", "r2"),
  resources = list(
    n_s1 = ~ clusters_s1 * m,
    n_s2 = ~ clusters_s2 * m,
    clusters_s1 = ~ clusters_s1,
    clusters_s2 = ~ clusters_s2,
    periods_s1 = ~ t1,
    periods_s2 = ~ t2
  ),
  cost_structure = list(
    weights = c(n = 1, clusters = "rho", periods = "rho"),
    stage2_resources = c("n_s2", "clusters_s2", "periods_s2")
  ),
  model_builder = staggered_model_builder,
  design_type = "staggered_enrollment"
)

# Step 3: Create design object
staggered_design <- structure(
  list(
    spec = staggered_spec,
    fixed_params = list(icc = 0.05, cac = 0.8, delta = 0.3),
    stage1_grid = data.frame(K = 30, t1 = 4, r1 = 4, m = 25),
    stage2_grid_fn = function(s1) expand.grid(t2 = 0:6, r2 = c(2, 4, 6, 8)),
    rho = 50
  ),
  class = c("staggered_crt", "crt_design")
)

# Step 4: Analyse
results <- adaptive_analysis(staggered_design, target_power = 0.8)
summary(results)
plot(results, type = "decision")
```

---

## Troubleshooting

### Common Issues

| Problem | Solution |
|---------|----------|
| "I1_eff is NA or Inf" | Check that idx1 and idx2 are correct; ensure X has full rank |
| "Model building is slow" | Implement caching with cache_key |
| "Costs seem wrong" | Check resource formulas match your model_builder outputs |
| "Sample sizes don't match" | Ensure n1, n2 in model_builder are actual observation counts |

### Debugging Tips

```r
# Test model builder directly
cache <- new.env()
test_output <- my_model_builder(
  design_params = data.frame(m2 = 30, k2 = 2),
  fixed_params = list(k1 = 10, m1 = 25, icc = 0.05, delta = 0.25, cache = cache)
)
print(test_output)

# Check efficient information is positive and sensible
stopifnot(test_output$I1_eff > 0)
stopifnot(test_output$I2_eff > 0)
stopifnot(abs(test_output$w1^2 + test_output$w2^2 - 1) < 0.001)
```
