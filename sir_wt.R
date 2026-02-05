#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Stochastic SIR renewal-process simulation with Wallinga-Teunis R_t estimation
#
# - R_t = 1.1 for t in [0, 50), R_t = 1.2 for t in [50, 100]
# - Generation interval: Gamma(mean=5, sd=2.5)
# - 200 simulations averaged
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)

# ── Parameters ───────────────────────────────────────────────────────────────
T_MAX      <- 100
N_SIMS     <- 200
GEN_MEAN   <- 5.0
GEN_SD     <- 2.5
INIT_CASES <- 10
MAX_LAG    <- 30

# True R_t schedule
true_rt <- function(t) ifelse(t < 50, 1.1, 1.2)

# Generation interval: discretised Gamma PMF
gen_shape <- (GEN_MEAN / GEN_SD)^2
gen_scale <- GEN_SD^2 / GEN_MEAN

gi_pmf <- sapply(1:MAX_LAG, function(s) {
  pgamma(s + 0.5, shape = gen_shape, scale = gen_scale) -
    pgamma(s - 0.5, shape = gen_shape, scale = gen_scale)
})
gi_pmf <- gi_pmf / sum(gi_pmf)

# ── Simulation: Renewal process ──────────────────────────────────────────────
simulate_renewal <- function(seed) {
  set.seed(seed)
  inc <- numeric(T_MAX + 1)
  inc[1] <- INIT_CASES          # index 1 = day 0

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

# ── Wallinga-Teunis estimator ────────────────────────────────────────────────
wallinga_teunis <- function(inc) {
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

      # denominator: all possible infectors of case at tj
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

# ── Run simulations ──────────────────────────────────────────────────────────
cat("Running 200 simulations...\n")

all_inc <- matrix(NA_real_, nrow = N_SIMS, ncol = T_MAX + 1)
all_rt  <- matrix(NA_real_, nrow = N_SIMS, ncol = T_MAX + 1)

successful <- 0
seed <- 0

while (successful < N_SIMS) {
  seed <- seed + 1
  inc  <- simulate_renewal(seed)

  # skip degenerate runs
  if (sum(inc[12:(T_MAX + 1)]) < 5) next

  successful <- successful + 1
  rt_est <- wallinga_teunis(inc)
  all_inc[successful, ] <- inc
  all_rt[successful, ]  <- rt_est

  if (successful %% 50 == 0) cat(sprintf("  %d/%d done\n", successful, N_SIMS))
}
cat("All simulations complete.\n")

# ── Summary statistics ───────────────────────────────────────────────────────
days <- 0:T_MAX

mean_inc <- colMeans(all_inc)
q05_inc  <- apply(all_inc, 2, quantile, 0.05)
q25_inc  <- apply(all_inc, 2, quantile, 0.25)
q75_inc  <- apply(all_inc, 2, quantile, 0.75)
q95_inc  <- apply(all_inc, 2, quantile, 0.95)

mean_rt <- colMeans(all_rt, na.rm = TRUE)
q05_rt  <- apply(all_rt, 2, quantile, 0.05, na.rm = TRUE)
q25_rt  <- apply(all_rt, 2, quantile, 0.25, na.rm = TRUE)
q75_rt  <- apply(all_rt, 2, quantile, 0.75, na.rm = TRUE)
q95_rt  <- apply(all_rt, 2, quantile, 0.95, na.rm = TRUE)

# Trim unreliable edges for WT panel
trim_start <- 4
trim_end   <- 6
idx <- (trim_start + 1):(T_MAX + 1 - trim_end)
d_trim <- days[idx]

# ── Plot ─────────────────────────────────────────────────────────────────────
png("/home/claude/sir_wt_plot_R.png", width = 10, height = 7.5, units = "in", res = 200)

par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 1), bg = "#FAFAF8", family = "serif",
    cex.lab = 1.1, cex.axis = 0.9)

# ── Panel A: Epidemic curve ──────────────────────────────────────────────────
plot(NULL, xlim = c(0, T_MAX), ylim = c(0, max(q95_inc) * 1.05),
     xlab = "", ylab = "Daily incidence", main = "", axes = FALSE)
title(main = expression(bold("A") ~ "  Stochastic SIR renewal process  (N = 200 simulations)"),
      adj = 0, cex.main = 1.25, line = 1)
axis(1); axis(2, las = 1)

# 90% band
polygon(c(days, rev(days)), c(q05_inc, rev(q95_inc)),
        col = adjustcolor("#BFD8E5", 0.5), border = NA)
# 50% band
polygon(c(days, rev(days)), c(q25_inc, rev(q75_inc)),
        col = adjustcolor("#7AB8D4", 0.6), border = NA)
# Mean
lines(days, mean_inc, col = "#1A5276", lwd = 2)

abline(v = 50, col = "#C0392B", lty = 2, lwd = 0.9)
text(52, max(mean_inc[41:61]) * 1.4, expression(italic(R[t]*": 1.1 -> 1.2")),
     col = "#C0392B", cex = 0.85, adj = 0)

# Legend - lines only
legend("topleft",
       legend = c("Mean incidence", "50% interval", "90% interval"),
       lwd = c(2, 10, 10),
       col = c("#1A5276", adjustcolor("#7AB8D4", 0.6), adjustcolor("#BFD8E5", 0.5)),
       bty = "n", cex = 0.8)

# ── Panel B: Wallinga-Teunis R_t ────────────────────────────────────────────
plot(NULL, xlim = c(0, T_MAX), ylim = c(0.5, 2.0),
     xlab = "Time (days)", ylab = expression(Reproduction ~ number ~ R[t]),
     main = "", axes = FALSE)
title(main = expression(bold("B") ~ "  Wallinga–Teunis reproduction number"),
      adj = 0, cex.main = 1.25, line = 1)
axis(1); axis(2, las = 1)

# 90% band
polygon(c(d_trim, rev(d_trim)), c(q05_rt[idx], rev(q95_rt[idx])),
        col = adjustcolor("#E8DAEF", 0.55), border = NA)
# 50% band
polygon(c(d_trim, rev(d_trim)), c(q25_rt[idx], rev(q75_rt[idx])),
        col = adjustcolor("#C39BD3", 0.6), border = NA)
# Mean WT
lines(d_trim, mean_rt[idx], col = "#6C3483", lwd = 2)

# True R_t (step function)
segments(0, 1.1, 50, 1.1, col = "#C0392B", lwd = 1.8)
segments(50, 1.2, 100, 1.2, col = "#C0392B", lwd = 1.8)

abline(h = 1.0, col = "#888888", lty = 3, lwd = 0.7)
abline(v = 50, col = "#C0392B", lty = 2, lwd = 0.9)

legend("topright",
       legend = c("Mean WT estimate", expression(True ~ R[t]),
                  "50% interval", "90% interval"),
       lwd = c(2, 1.8, 10, 10),
       col = c("#6C3483", "#C0392B",
               adjustcolor("#C39BD3", 0.6), adjustcolor("#E8DAEF", 0.55)),
       bty = "n", cex = 0.8)

# Footer annotation
mtext(paste0("Generation interval: Gamma(mu=", GEN_MEAN, ", sd=", GEN_SD,
             ")  |  Initial cases: ", INIT_CASES, "  |  Simulations: ", N_SIMS),
      side = 1, outer = TRUE, line = -1.2, cex = 0.75, col = "#666666", font = 3)

dev.off()
cat("Plot saved.\n")
