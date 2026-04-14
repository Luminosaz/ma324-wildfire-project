## ============================================================
## test_samplers.R
## Validation tests for step3_samplers.R
## ============================================================

source("step3_samplers.R")
set.seed(1)

cat("============================================================\n")
cat("  TEST 1: dwind_dir — density validation\n")
cat("============================================================\n")

## density must integrate to 1
integral = integrate(dwind_dir, lower = 0, upper = 2 * pi)$value
cat("Integral of f over [0, 2pi):", integral, "\n")
stopifnot(abs(integral - 1) < 1e-6)

## density must be non-negative everywhere
theta_grid = seq(0, 2 * pi, length.out = 10000)
f_vals     = dwind_dir(theta_grid)
cat("Min density value:          ", min(f_vals), "\n")
stopifnot(all(f_vals >= 0))

cat("PASSED\n\n")


cat("============================================================\n")
cat("  TEST 2: sample_wind_dir — acceptance-rejection sampler\n")
cat("============================================================\n")

n_dir = 10000
res   = sample_wind_dir(n_dir)

cat("n requested:      ", n_dir, "\n")
cat("n returned:       ", length(res$samples), "\n")
cat("M (envelope):     ", round(res$M, 6), "\n")
cat("Acceptance rate:  ", round(res$acceptance_rate, 4), "\n")
cat("Theoretical rate: ", round(res$theoretical_rate, 4), "\n")

## correct number of samples returned
stopifnot(length(res$samples) == n_dir)

## all samples in [0, 2pi)
stopifnot(all(res$samples >= 0 & res$samples < 2 * pi))

## theoretical_rate matches formula
stopifnot(abs(res$theoretical_rate - 1 / (2 * pi * res$M)) < 1e-10)

## acceptance rate within reasonable range of theoretical
cat("Rate ratio (actual/theoretical):",
    round(res$acceptance_rate / res$theoretical_rate, 4), "\n")

## KS test: compare samples against mixture CDF (via numerical integration)
cdf_vals = sapply(sort(res$samples), function(x) {
  integrate(dwind_dir, lower = 0, upper = x)$value
})
ks = ks.test(cdf_vals, "punif")
cat("KS test p-value:  ", round(ks$p.value, 4),
    "(fail to reject if > 0.05)\n")
stopifnot(ks$p.value > 0.01)

cat("PASSED\n\n")


cat("============================================================\n")
cat("  TEST 3: sample_wind_speed — Weibull(shape=2, scale=7)\n")
cat("============================================================\n")

n_ws = 10000
ws   = sample_wind_speed(n_ws)

## theoretical moments for Weibull(k=2, A=7)
theo_mean = 7 * gamma(1 + 1 / 2)
theo_var  = 7^2 * (gamma(1 + 2 / 2) - gamma(1 + 1 / 2)^2)

cat("n returned:       ", length(ws), "\n")
cat("Sample mean:      ", round(mean(ws), 3),
    " (theoretical:", round(theo_mean, 3), ")\n")
cat("Sample var:       ", round(var(ws), 3),
    " (theoretical:", round(theo_var, 3), ")\n")

## correct length
stopifnot(length(ws) == n_ws)

## all non-negative (Weibull support)
stopifnot(all(ws >= 0))

## KS test against true Weibull CDF
ks_ws = ks.test(ws, "pweibull", shape = 2, scale = 7)
cat("KS test p-value:  ", round(ks_ws$p.value, 4), "\n")
stopifnot(ks_ws$p.value > 0.01)

cat("PASSED\n\n")


cat("============================================================\n")
cat("  TEST 4: sample_ignition — uniform over vegetated cells\n")
cat("============================================================\n")

landscape = as.matrix(read.table(
  "../../Project data and resources/landscape.csv",
  sep = ",", header = FALSE))
cat("Landscape dims:   ", nrow(landscape), "x", ncol(landscape), "\n")

n_veg   = sum(landscape != 0)
n_bare  = sum(landscape == 0)
cat("Vegetated cells:  ", n_veg, "\n")
cat("Bare cells:       ", n_bare, "\n")

n_ign = 10000
ign   = sample_ignition(landscape, n_ign)

cat("n returned:       ", nrow(ign), "\n")
cat("Columns:          ", paste(colnames(ign), collapse = ", "), "\n")

## correct dimensions and column names
stopifnot(nrow(ign) == n_ign)
stopifnot(ncol(ign) == 2)
stopifnot(all(colnames(ign) == c("row", "col")))

## all indices within grid bounds
stopifnot(all(ign[, "row"] >= 1 & ign[, "row"] <= nrow(landscape)))
stopifnot(all(ign[, "col"] >= 1 & ign[, "col"] <= ncol(landscape)))

## every sampled cell is vegetated (nonzero)
cell_vals = landscape[ign]
cat("All vegetated?    ", all(cell_vals != 0), "\n")
stopifnot(all(cell_vals != 0))

## chi-squared test for uniformity over vegetated cells
veg_idx    = which(landscape != 0, arr.ind = TRUE)
cell_id    = function(rc) (rc[, 1] - 1) * ncol(landscape) + rc[, 2]
sampled_id = cell_id(ign)
freq       = table(factor(sampled_id, levels = cell_id(veg_idx)))
chi2       = chisq.test(as.numeric(freq))
cat("Chi-sq p-value:   ", round(chi2$p.value, 4),
    "(fail to reject uniformity if > 0.05)\n")
stopifnot(chi2$p.value > 0.01)

cat("PASSED\n\n")


cat("============================================================\n")
cat("  ALL TESTS PASSED\n")
cat("============================================================\n")
