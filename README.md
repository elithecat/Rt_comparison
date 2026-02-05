Wallinga-Teunis R<sub>t</sub> Estimation: A Simulation Study
================
Epidemiological Modeling
2026-02-05

-   [Introduction](#introduction)
    -   [The Reproduction Number](#the-reproduction-number)
    -   [The Wallinga-Teunis Method](#the-wallinga-teunis-method)
-   [Simulation Setup](#simulation-setup)
    -   [Generation Interval
        Distribution](#generation-interval-distribution)
    -   [True $R_t$ Schedule](#true-r_t-schedule)
-   [Renewal Process Simulation](#renewal-process-simulation)
-   [Custom Wallinga-Teunis
    Implementation](#custom-wallinga-teunis-implementation)
-   [EpiEstim Wallinga-Teunis
    Implementation](#epiestim-wallinga-teunis-implementation)
-   [Running the Simulations](#running-the-simulations)
-   [Results](#results)
    -   [Summary Statistics](#summary-statistics)
    -   [Epidemic Curves](#epidemic-curves)
    -   [Comparison of $R_t$ Estimators](#comparison-of-r_t-estimators)
    -   [Direct Method Comparison](#direct-method-comparison)
    -   [Numerical Comparison](#numerical-comparison)
-   [Discussion](#discussion)
    -   [Key Findings](#key-findings)
    -   [Wallinga-Teunis Properties](#wallinga-teunis-properties)
    -   [When to Use Each
        Implementation](#when-to-use-each-implementation)
-   [References](#references)
-   [Session Information](#session-information)

## Introduction

This document demonstrates the estimation of the time-varying
reproduction number $R_t$ using the **Wallinga-Teunis method**. We
compare a custom implementation with the `EpiEstim` package
implementation using simulated epidemic data.

### The Reproduction Number

The reproduction number $R_t$ represents the average number of secondary
infections caused by a single infected individual at time $t$. It is a
key metric for understanding epidemic dynamics:

-   $R_t > 1$: Epidemic is growing
-   $R_t = 1$: Epidemic is stable
-   $R_t < 1$: Epidemic is declining

### The Wallinga-Teunis Method

The Wallinga-Teunis method (Wallinga & Teunis, 2004) estimates the
**case reproduction number** by:

1.  Computing the relative likelihood that case $i$ infected case $j$
    based on the serial interval distribution
2.  Summing these likelihoods for each potential infector to estimate
    their individual reproduction number
3.  Averaging across cases occurring at the same time

The key formula is:

$$R_i = \sum_{j} p_{ij}$$

where $p_{ij}$ is the probability that case $i$ infected case $j$,
computed as:

$$p_{ij} = \frac{w(t_j - t_i)}{\sum_k w(t_j - t_k)}$$

and $w(\cdot)$ is the serial interval distribution.

## Simulation Setup

We simulate epidemics using a **renewal process** with known $R_t$
values:

| Parameter           | Value                        | Description                         |
|---------------------|------------------------------|-------------------------------------|
| $R_t$ (days 0-49)   | 1.1                          | Reproduction number in first phase  |
| $R_t$ (days 50-100) | 1.2                          | Reproduction number in second phase |
| Generation interval | Gamma($\mu$=5, $\sigma$=2.5) | Time between infections             |
| Initial cases       | 10                           | Seed infections at $t=0$            |
| Simulations         | 200                          | Number of stochastic realizations   |

``` r
# Simulation parameters
T_MAX      <- 100    # Total simulation days
N_SIMS     <- 200    # Number of simulations
GEN_MEAN   <- 5.0    # Generation interval mean (days)
GEN_SD     <- 2.5    # Generation interval SD (days)
INIT_CASES <- 10     # Initial number of cases
MAX_LAG    <- 30     # Maximum generation interval lag

set.seed(42)
```

### Generation Interval Distribution

The generation interval (time from infection of a primary case to
infection of a secondary case) follows a Gamma distribution, which we
discretize for the renewal process.

``` r
# Gamma distribution parameters
gen_shape <- (GEN_MEAN / GEN_SD)^2
gen_scale <- GEN_SD^2 / GEN_MEAN

# Discretized PMF
gi_pmf <- sapply(1:MAX_LAG, function(s) {
  pgamma(s + 0.5, shape = gen_shape, scale = gen_scale) -
    pgamma(s - 0.5, shape = gen_shape, scale = gen_scale)
})
gi_pmf <- gi_pmf / sum(gi_pmf)

# Plot the generation interval distribution
par(mar = c(4.5, 4.5, 2, 1), bg = "#FAFAF8")
barplot(gi_pmf, names.arg = 1:MAX_LAG, col = "#3498DB", border = NA,
        xlab = "Days since infection", ylab = "Probability",
        main = "Discretized Generation Interval Distribution")
abline(h = 0, col = "gray50")
text(GEN_MEAN * 1.2, max(gi_pmf) * 0.9,
     bquote(mu == .(GEN_MEAN) ~ "days," ~ sigma == .(GEN_SD) ~ "days"),
     adj = 0, cex = 0.9)
```

<img src="/workspace/figures/generation-interval-1.png" alt="" style="display: block; margin: auto;" />

### True $R_t$ Schedule

``` r
# True R_t as a function of time
true_rt <- function(t) ifelse(t < 50, 1.1, 1.2)

# Visualize the true R_t
par(mar = c(4.5, 4.5, 2, 1), bg = "#FAFAF8")
curve(sapply(x, true_rt), from = 0, to = 100, n = 1000,
      xlab = "Time (days)", ylab = expression(R[t]),
      main = expression("True " * R[t] * " Schedule"),
      col = "#C0392B", lwd = 2, ylim = c(0.8, 1.5))
abline(h = 1, lty = 3, col = "gray50")
abline(v = 50, lty = 2, col = "#C0392B")
```

<img src="/workspace/figures/true-rt-1.png" alt="" style="display: block; margin: auto;" />

## Renewal Process Simulation

The renewal equation generates new infections based on past incidence
and the generation interval:

$$I_t = \text{Poisson}\left(\sum_{s=1}^{t} I_{t-s} \cdot R_{t-s} \cdot w_s\right)$$

``` r
simulate_renewal <- function(seed) {
  set.seed(seed)
  inc <- numeric(T_MAX + 1)
  inc[1] <- INIT_CASES  # index 1 = day 0

  for (t in 1:T_MAX) {
    lam <- 0
    for (lag in 1:min(t, MAX_LAG)) {
      rt_val <- true_rt(t - lag)
      lam <- lam + inc[t - lag + 1] * rt_val * gi_pmf[lag]
    }
    inc[t + 1] <- rpois(1, max(lam, 0))
  }
  inc
}

# Example single simulation
example_inc <- simulate_renewal(1)
par(mar = c(4.5, 4.5, 2, 1), bg = "#FAFAF8")
barplot(example_inc, names.arg = 0:T_MAX, col = "#2ECC71", border = NA,
        xlab = "Time (days)", ylab = "Daily incidence",
        main = "Example Epidemic Curve (Single Simulation)")
abline(v = 50 * 1.2, lty = 2, col = "#C0392B")
```

<img src="/workspace/figures/simulate-renewal-1.png" alt="" style="display: block; margin: auto;" />

## Custom Wallinga-Teunis Implementation

The custom implementation directly computes the case reproduction number
following the original method.

``` r
wallinga_teunis_custom <- function(inc) {
  TT <- length(inc)
  Rt <- rep(NA_real_, TT)

  for (ti in 1:TT) {
    if (inc[ti] == 0) next
    total <- 0

    tj_max <- min(ti + MAX_LAG, TT)
    if (ti + 1 > tj_max) next

    for (tj in (ti + 1):tj_max) {
      lag_ij <- tj - ti
      w_ij   <- gi_pmf[lag_ij]

      # Denominator: sum over all possible infectors of case at tj
      denom <- 0
      for (tk in max(1, tj - MAX_LAG):(tj - 1)) {
        lag_kj <- tj - tk
        if (lag_kj >= 1 && lag_kj <= MAX_LAG) {
          denom <- denom + inc[tk] * gi_pmf[lag_kj]
        }
      }

      if (denom > 0) {
        total <- total + inc[tj] * w_ij / denom
      }
    }
    Rt[ti] <- total
  }
  Rt
}
```

## EpiEstim Wallinga-Teunis Implementation

The `EpiEstim` package provides a well-tested implementation of the
Wallinga-Teunis method with additional features for uncertainty
quantification.

``` r
library(EpiEstim)

estimate_wt_epiestim <- function(inc) {
  # Serial interval distribution for EpiEstim (starts at lag 0)
  si_distr <- c(0, gi_pmf)

  # Define time windows (single-day windows)
  n_days <- length(inc)
  t_start <- 2:n_days
  t_end <- t_start

  # Run EpiEstim's wallinga_teunis
  wt_result <- wallinga_teunis(
    incid = inc,
    method = "non_parametric_si",
    config = list(
      t_start = t_start,
      t_end = t_end,
      si_distr = si_distr,
      n_sim = 100
    )
  )

  # Extract mean R estimates

rt_est <- rep(NA_real_, length(inc))
  if (!is.null(wt_result$R)) {
    for (i in seq_len(nrow(wt_result$R))) {
      t_idx <- wt_result$R$t_start[i]
      if (t_idx >= 1 && t_idx <= length(inc)) {
        rt_est[t_idx] <- wt_result$R$`Mean(R)`[i]
      }
    }
  }
  rt_est
}
```

## Running the Simulations

We run 200 simulations and estimate $R_t$ using both methods.

``` r
# Storage matrices
all_inc       <- matrix(NA_real_, nrow = N_SIMS, ncol = T_MAX + 1)
all_rt_custom <- matrix(NA_real_, nrow = N_SIMS, ncol = T_MAX + 1)
all_rt_epi    <- matrix(NA_real_, nrow = N_SIMS, ncol = T_MAX + 1)

successful <- 0
seed <- 0

cat("Running simulations...\n")
```

    ## Running simulations...

``` r
while (successful < N_SIMS) {
  seed <- seed + 1
  inc <- simulate_renewal(seed)

  # Skip degenerate runs (near-extinction)
  if (sum(inc[12:(T_MAX + 1)]) < 5) next

  successful <- successful + 1

  # Estimate R_t with both methods
  rt_custom <- wallinga_teunis_custom(inc)
  rt_epi <- tryCatch(
    estimate_wt_epiestim(inc),
    error = function(e) rep(NA_real_, length(inc))
  )

  all_inc[successful, ]       <- inc
  all_rt_custom[successful, ] <- rt_custom
  all_rt_epi[successful, ]    <- rt_epi

  if (successful %% 50 == 0) cat(sprintf("  %d/%d complete\n", successful, N_SIMS))
}
```

    ##   50/200 complete
    ##   100/200 complete
    ##   150/200 complete
    ##   200/200 complete

``` r
cat("All simulations complete.\n")
```

    ## All simulations complete.

## Results

### Summary Statistics

``` r
days <- 0:T_MAX

# Incidence statistics
mean_inc <- colMeans(all_inc)
q05_inc  <- apply(all_inc, 2, quantile, 0.05)
q25_inc  <- apply(all_inc, 2, quantile, 0.25)
q75_inc  <- apply(all_inc, 2, quantile, 0.75)
q95_inc  <- apply(all_inc, 2, quantile, 0.95)

# Custom WT statistics
mean_rt_custom <- colMeans(all_rt_custom, na.rm = TRUE)
q05_rt_custom  <- apply(all_rt_custom, 2, quantile, 0.05, na.rm = TRUE)
q25_rt_custom  <- apply(all_rt_custom, 2, quantile, 0.25, na.rm = TRUE)
q75_rt_custom  <- apply(all_rt_custom, 2, quantile, 0.75, na.rm = TRUE)
q95_rt_custom  <- apply(all_rt_custom, 2, quantile, 0.95, na.rm = TRUE)

# EpiEstim WT statistics
mean_rt_epi <- colMeans(all_rt_epi, na.rm = TRUE)
q05_rt_epi  <- apply(all_rt_epi, 2, quantile, 0.05, na.rm = TRUE)
q25_rt_epi  <- apply(all_rt_epi, 2, quantile, 0.25, na.rm = TRUE)
q75_rt_epi  <- apply(all_rt_epi, 2, quantile, 0.75, na.rm = TRUE)
q95_rt_epi  <- apply(all_rt_epi, 2, quantile, 0.95, na.rm = TRUE)

# Trim edges where estimates are unreliable
trim_start <- 4
trim_end   <- 6
idx <- (trim_start + 1):(T_MAX + 1 - trim_end)
d_trim <- days[idx]
```

### Epidemic Curves

``` r
par(mar = c(4.5, 4.5, 3, 1), bg = "#FAFAF8")

plot(NULL, xlim = c(0, T_MAX), ylim = c(0, max(q95_inc) * 1.05),
     xlab = "Time (days)", ylab = "Daily incidence", axes = FALSE)
title(main = "Stochastic SIR Renewal Process (N = 200 simulations)",
      cex.main = 1.2)
axis(1); axis(2, las = 1)

# 90% interval
polygon(c(days, rev(days)), c(q05_inc, rev(q95_inc)),
        col = adjustcolor("#BFD8E5", 0.5), border = NA)
# 50% interval
polygon(c(days, rev(days)), c(q25_inc, rev(q75_inc)),
        col = adjustcolor("#7AB8D4", 0.6), border = NA)
# Mean
lines(days, mean_inc, col = "#1A5276", lwd = 2)

abline(v = 50, col = "#C0392B", lty = 2)
text(52, max(mean_inc) * 0.9, expression(R[t] * ": 1.1 -> 1.2"),
     col = "#C0392B", cex = 0.85, adj = 0)

legend("topleft",
       legend = c("Mean incidence", "50% interval", "90% interval"),
       lwd = c(2, 10, 10),
       col = c("#1A5276", adjustcolor("#7AB8D4", 0.6), adjustcolor("#BFD8E5", 0.5)),
       bty = "n", cex = 0.85)
```

<img src="/workspace/figures/plot-incidence-1.png" alt="" style="display: block; margin: auto;" />

### Comparison of $R_t$ Estimators

``` r
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 1), bg = "#FAFAF8")

# Panel A: Custom Implementation
plot(NULL, xlim = c(0, T_MAX), ylim = c(0.5, 2.0),
     xlab = "Time (days)", ylab = expression(R[t]), axes = FALSE)
title(main = "Custom Wallinga-Teunis Implementation", cex.main = 1.2)
axis(1); axis(2, las = 1)

polygon(c(d_trim, rev(d_trim)), c(q05_rt_custom[idx], rev(q95_rt_custom[idx])),
        col = adjustcolor("#E8DAEF", 0.55), border = NA)
polygon(c(d_trim, rev(d_trim)), c(q25_rt_custom[idx], rev(q75_rt_custom[idx])),
        col = adjustcolor("#C39BD3", 0.6), border = NA)
lines(d_trim, mean_rt_custom[idx], col = "#6C3483", lwd = 2)

segments(0, 1.1, 50, 1.1, col = "#C0392B", lwd = 1.8)
segments(50, 1.2, 100, 1.2, col = "#C0392B", lwd = 1.8)
abline(h = 1.0, col = "#888888", lty = 3)
abline(v = 50, col = "#C0392B", lty = 2)

legend("topright",
       legend = c("Mean estimate", expression(True ~ R[t]), "50% interval", "90% interval"),
       lwd = c(2, 1.8, 10, 10),
       col = c("#6C3483", "#C0392B", adjustcolor("#C39BD3", 0.6), adjustcolor("#E8DAEF", 0.55)),
       bty = "n", cex = 0.8)

# Panel B: EpiEstim Implementation
plot(NULL, xlim = c(0, T_MAX), ylim = c(0.5, 2.0),
     xlab = "Time (days)", ylab = expression(R[t]), axes = FALSE)
title(main = "EpiEstim Wallinga-Teunis Implementation", cex.main = 1.2)
axis(1); axis(2, las = 1)

polygon(c(d_trim, rev(d_trim)), c(q05_rt_epi[idx], rev(q95_rt_epi[idx])),
        col = adjustcolor("#D5F5E3", 0.55), border = NA)
polygon(c(d_trim, rev(d_trim)), c(q25_rt_epi[idx], rev(q75_rt_epi[idx])),
        col = adjustcolor("#82E0AA", 0.6), border = NA)
lines(d_trim, mean_rt_epi[idx], col = "#1E8449", lwd = 2)

segments(0, 1.1, 50, 1.1, col = "#C0392B", lwd = 1.8)
segments(50, 1.2, 100, 1.2, col = "#C0392B", lwd = 1.8)
abline(h = 1.0, col = "#888888", lty = 3)
abline(v = 50, col = "#C0392B", lty = 2)

legend("topright",
       legend = c("Mean estimate", expression(True ~ R[t]), "50% interval", "90% interval"),
       lwd = c(2, 1.8, 10, 10),
       col = c("#1E8449", "#C0392B", adjustcolor("#82E0AA", 0.6), adjustcolor("#D5F5E3", 0.55)),
       bty = "n", cex = 0.8)
```

<img src="/workspace/figures/plot-comparison-1.png" alt="" style="display: block; margin: auto;" />

### Direct Method Comparison

``` r
par(mar = c(4.5, 4.5, 3, 1), bg = "#FAFAF8")

plot(NULL, xlim = c(0, T_MAX), ylim = c(0.8, 1.6),
     xlab = "Time (days)", ylab = expression(R[t]), axes = FALSE)
title(main = "Comparison: Custom vs EpiEstim Wallinga-Teunis", cex.main = 1.2)
axis(1); axis(2, las = 1)

# True R_t
segments(0, 1.1, 50, 1.1, col = "#C0392B", lwd = 2)
segments(50, 1.2, 100, 1.2, col = "#C0392B", lwd = 2)

# Custom estimate
lines(d_trim, mean_rt_custom[idx], col = "#6C3483", lwd = 2)

# EpiEstim estimate
lines(d_trim, mean_rt_epi[idx], col = "#1E8449", lwd = 2, lty = 2)

abline(h = 1.0, col = "#888888", lty = 3)
abline(v = 50, col = "#C0392B", lty = 2, lwd = 0.8)

legend("topright",
       legend = c(expression(True ~ R[t]), "Custom WT", "EpiEstim WT"),
       lwd = c(2, 2, 2),
       lty = c(1, 1, 2),
       col = c("#C0392B", "#6C3483", "#1E8449"),
       bty = "n", cex = 0.9)
```

<img src="/workspace/figures/plot-overlay-1.png" alt="" style="display: block; margin: auto;" />

### Numerical Comparison

``` r
# Calculate mean absolute error for each method
mae_custom <- mean(abs(mean_rt_custom[idx] - sapply(d_trim, true_rt)), na.rm = TRUE)
mae_epi    <- mean(abs(mean_rt_epi[idx] - sapply(d_trim, true_rt)), na.rm = TRUE)

# Correlation between methods
cor_methods <- cor(mean_rt_custom[idx], mean_rt_epi[idx], use = "complete.obs")

comparison_table <- data.frame(
  Metric = c("Mean Absolute Error", "Correlation with True R_t", "Method Correlation"),
  Custom_WT = c(round(mae_custom, 4),
                round(cor(mean_rt_custom[idx], sapply(d_trim, true_rt), use = "complete.obs"), 4),
                "-"),
  EpiEstim_WT = c(round(mae_epi, 4),
                  round(cor(mean_rt_epi[idx], sapply(d_trim, true_rt), use = "complete.obs"), 4),
                  round(cor_methods, 4))
)

knitr::kable(comparison_table,
             col.names = c("Metric", "Custom WT", "EpiEstim WT"),
             caption = "Performance Comparison of Wallinga-Teunis Implementations",
             align = "lcc")
```

| Metric                    | Custom WT | EpiEstim WT |
|:--------------------------|:---------:|:-----------:|
| Mean Absolute Error       |  0.0219   |   0.0219    |
| Correlation with True R_t |  0.6798   |   0.6798    |
| Method Correlation        |    \-     |   1.0000    |

Performance Comparison of Wallinga-Teunis Implementations

## Discussion

### Key Findings

1.  **Both implementations recover the true $R_t$**: The mean estimates
    closely track the true values of 1.1 (days 0-49) and 1.2 (days
    50-100).

2.  **Methods produce nearly identical results**: The high correlation
    between custom and EpiEstim implementations validates the custom
    code.

3.  **Characteristic lag at change point**: Both methods show a brief
    transition period around day 50, which is expected behavior for the
    Wallinga-Teunis estimator.

4.  **Edge effects**: Estimates are less reliable at the beginning and
    end of the time series, hence the trimming in visualizations.

### Wallinga-Teunis Properties

The Wallinga-Teunis method estimates the **case reproduction number**,
which:

-   Assigns $R$ to the time of the *infector* (not the infectee)
-   Is a retrospective measure (requires future data)
-   Has lower variance than instantaneous $R_t$ estimators
-   Shows a characteristic rightward shift relative to the “true” $R_t$

### When to Use Each Implementation

| Scenario                   | Recommendation                                    |
|----------------------------|---------------------------------------------------|
| Production analysis        | Use `EpiEstim` (well-tested, maintained)          |
| Custom modifications       | Start with custom implementation                  |
| Teaching/learning          | Custom implementation shows the algorithm clearly |
| Uncertainty quantification | `EpiEstim` provides built-in CI methods           |

## References

1.  Wallinga, J., & Teunis, P. (2004). Different epidemic curves for
    severe acute respiratory syndrome reveal similar impacts of control
    measures. *American Journal of Epidemiology*, 160(6), 509-516.

2.  Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013). A new
    framework and software to estimate time-varying reproduction numbers
    during epidemics. *American Journal of Epidemiology*, 178(9),
    1505-1512.

3.  EpiEstim R package: <https://cran.r-project.org/package=EpiEstim>

## Session Information

``` r
sessionInfo()
```

    ## R version 4.2.2 Patched (2022-11-10 r83330)
    ## Platform: aarch64-unknown-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 12 (bookworm)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/aarch64-linux-gnu/blas/libblas.so.3.11.0
    ## LAPACK: /usr/lib/aarch64-linux-gnu/lapack/liblapack.so.3.11.0
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] EpiEstim_2.2-5
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.1.1            pillar_1.11.1         RColorBrewer_1.1-3   
    ##  [4] compiler_4.2.2        plyr_1.8.9            tools_4.2.2          
    ##  [7] fitdistrplus_1.2-6    digest_0.6.39         tibble_3.3.1         
    ## [10] evaluate_1.0.5        lifecycle_1.0.5       gtable_0.3.6         
    ## [13] lattice_0.20-45       pkgconfig_2.0.3       rlang_1.1.7          
    ## [16] Matrix_1.5-4          cli_3.6.5             yaml_2.3.12          
    ## [19] SparseM_1.84-2        xfun_0.56             fastmap_1.2.0        
    ## [22] gridExtra_2.3         coda_0.19-4.1         dplyr_1.1.4          
    ## [25] stringr_1.6.0         knitr_1.51            generics_0.1.4       
    ## [28] vctrs_0.7.1           coarseDataTools_0.7.2 MatrixModels_0.5-1   
    ## [31] tidyselect_1.2.1      grid_4.2.2            incidence_1.7.6      
    ## [34] glue_1.8.0            R6_2.6.1              otel_0.2.0           
    ## [37] survival_3.8-6        rmarkdown_2.30        farver_2.1.2         
    ## [40] reshape2_1.4.5        ggplot2_4.0.1         magrittr_2.0.4       
    ## [43] codetools_0.2-19      scales_1.4.0          htmltools_0.5.9      
    ## [46] mcmc_0.9-8            MASS_7.3-58.2         splines_4.2.2        
    ## [49] S7_0.2.1              quantreg_5.95         stringi_1.8.7        
    ## [52] MCMCpack_1.7-1
