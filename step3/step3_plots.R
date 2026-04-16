## ============================================================
## step3_plots.R
## Three multi-panel figures for Step 3 report.
## ============================================================

source("step3_samplers.R")   # for dwind_dir

scenarios = read.csv("step3_scenarios_N5000.csv")
risk      = read.csv("step3_risk_measures_N5000.csv")

fig_dir = "../figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## colour palette (consistent across figures)
plan_cols  = c(baseline = "#E64B35", corridor = "#4DBBD5", proposed = "#00A087")
plan_names = c(baseline = "No firebreaks", corridor = "Corridor", proposed = "Proposed (B=13)")


## ============================================================
## FIGURE 1: Sampler validation (2 x 1)
## Uses the same N = 5000 batch consumed by the Monte Carlo so that
## the acceptance / KS numbers in the Results match the figure.
## ============================================================
set.seed(1)
wd_samples = sample_wind_dir(5000)$samples
ws_samples = sample_wind_speed(5000)

png(file.path(fig_dir, "fig1_sampler_validation.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), cex.lab = 1.1)

## — panel A: wind direction —
hist(wd_samples, breaks = 60, freq = FALSE, col = "grey85", border = "white",
     main = "Wind direction sampler",
     xlab = expression(theta ~ "(radians)"), ylab = "Density",
     xlim = c(0, 2 * pi))
theta_seq = seq(0, 2 * pi, length.out = 500)
lines(theta_seq, dwind_dir(theta_seq), col = "#E64B35", lwd = 2.5)
legend("topright", legend = c("Samples", "True density"),
       fill = c("grey85", NA), border = c("grey60", NA),
       lty = c(NA, 1), lwd = c(NA, 2.5), col = c(NA, "#E64B35"),
       bty = "n", cex = 0.9)

## — panel B: wind speed —
hist(ws_samples, breaks = 50, freq = FALSE, col = "grey85", border = "white",
     main = "Wind speed sampler",
     xlab = "Wind speed (km/h)", ylab = "Density")
v_seq = seq(0, max(ws_samples) * 1.05, length.out = 500)
lines(v_seq, dweibull(v_seq, shape = 2, scale = 7), col = "#4DBBD5", lwd = 2.5)
legend("topright", legend = c("Samples", "Weibull(2, 7)"),
       fill = c("grey85", NA), border = c("grey60", NA),
       lty = c(NA, 1), lwd = c(NA, 2.5), col = c(NA, "#4DBBD5"),
       bty = "n", cex = 0.9)

dev.off()
cat("Saved fig1_sampler_validation.png\n")


## ============================================================
## FIGURE 2: Damage densities by plan (1 x 3)
## ============================================================
x_max = quantile(scenarios$damage, 0.995)

png(file.path(fig_dir, "fig2_damage_densities.png"),
    width = 12, height = 4.5, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1), cex.lab = 1.1)

for (pname in c("baseline", "corridor", "proposed")) {
  dmg  = scenarios$damage[scenarios$plan == pname]
  rmsk = risk[risk$plan == pname, ]
  mu   = rmsk$estimate[rmsk$measure == "mean_damage"]
  cvar = rmsk$estimate[rmsk$measure == "CVaR_0.9"]
  col  = plan_cols[pname]

  hist(dmg, breaks = 50, freq = FALSE, col = adjustcolor(col, 0.3),
       border = "white", main = plan_names[pname],
       xlab = "Damage", ylab = "Density", xlim = c(0, x_max))

  ## mean + CVaR lines
  abline(v = mu,   col = "grey30", lwd = 2, lty = 2)
  abline(v = cvar, col = col,      lwd = 2, lty = 1)

  legend("topright",
         legend = c(sprintf("Mean = %.1f", mu),
                    sprintf("CVaR(0.9) = %.1f", cvar)),
         lty = c(2, 1), lwd = 2, col = c("grey30", col),
         bty = "n", cex = 0.85)
}

dev.off()
cat("Saved fig2_damage_densities.png\n")


## ============================================================
## FIGURE 3: Proposed plan diagnostics (1 x 2)
## ============================================================
prop = scenarios[scenarios$plan == "proposed", ]
cvar_threshold = quantile(prop$damage, 0.9)
is_worst = prop$damage >= cvar_threshold

png(file.path(fig_dir, "fig3_proposed_diagnostics.png"),
    width = 11, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1.5), cex.lab = 1.1)

## — panel A: damage vs wind speed, coloured by compass direction —
## bin wind_dir into 4 compass quadrants (wind comes FROM)
compass = cut(prop$wind_dir %% (2 * pi),
              breaks = c(0, pi/4, 3*pi/4, 5*pi/4, 7*pi/4, 2*pi),
              labels = c("N", "E", "S", "W", "N"),
              include.lowest = TRUE)
compass_cols = c(N = "#999999", E = "#E64B35", S = "#4DBBD5", W = "#00A087")
pt_col = adjustcolor(compass_cols[as.character(compass)], 0.5)

## only plot damage > 0 so the pattern is visible
show = prop$damage > 0
plot(prop$wind_speed[show & !is_worst], prop$damage[show & !is_worst],
     pch = 16, cex = 0.6, col = pt_col[show & !is_worst],
     xlim = range(prop$wind_speed), ylim = c(0, max(prop$damage)),
     xlab = "Wind speed (km/h)", ylab = "Damage",
     main = "Damage vs wind speed (proposed, damage > 0)")
points(prop$wind_speed[is_worst], prop$damage[is_worst],
       pch = 1, cex = 0.8, lwd = 1.2, col = "black")
abline(h = cvar_threshold, lty = 2, col = "grey40", lwd = 1.5)

legend("topleft",
       legend = c("N wind", "E wind", "S wind", "W wind", "Worst 10%"),
       pch = c(16, 16, 16, 16, 1),
       col = c(compass_cols, "black"),
       bty = "n", cex = 0.8, pt.cex = 1)

## — panel B: mean damage by ignition location (heatmap) —
agg = aggregate(damage ~ ign_row + ign_col, data = prop, FUN = mean)

dmg_grid = matrix(NA, nrow = 21, ncol = 21)
for (k in seq_len(nrow(agg))) {
  dmg_grid[agg$ign_row[k], agg$ign_col[k]] = agg$damage[k]
}

heat_cols = colorRampPalette(c("#FFFFCC", "#FD8D3C", "#BD0026"))(100)

image(x = 1:21, y = 1:21, z = t(dmg_grid[21:1, ]),
      col = heat_cols, useRaster = TRUE,
      xlab = "Column", ylab = "Row (1 = north)",
      main = "Mean damage by ignition cell")

rect(3 - 0.5,  21 - 6 + 0.5,  6 + 0.5,  21 - 3 + 1.5,
     border = "#00A087", lwd = 2)   # wetland NW
rect(16 - 0.5, 21 - 6 + 0.5, 19 + 0.5, 21 - 3 + 1.5,
     border = "#4DBBD5", lwd = 2)   # settlement NE

legend("bottomright",
       legend = c("Wetland", "Settlement"),
       border = c("#00A087", "#4DBBD5"), fill = NA,
       bty = "n", cex = 0.85)

dev.off()
cat("Saved fig3_proposed_diagnostics.png\n")
