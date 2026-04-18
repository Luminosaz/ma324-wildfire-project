# ============================================================
# run_R4_N200.R — extra sensitivity: (B=13, R=4) at N=200
# ============================================================
# Re-runs the main (B=13, R=4) configuration at N=200 so that
# the R sensitivity panel in Results compares all R values at
# a common N. Saves to sensitivity_R4_N200.rds (does NOT
# touch final_N1000.rds or any existing sensitivity_*.rds).
# ============================================================

source("../.Rprofile")

data_dir = "../../Project data and resources"
source(file.path(data_dir, "fire-simulator.r"))
source("../step3/step3_samplers.R")
source("../step3/step3_monte_carlo.R")
source("evaluate_all_plans.R")

landscape = as.matrix(read.csv(file.path(data_dir, "landscape.csv"),
                               header = FALSE))
targets = read.csv(file.path(data_dir, "targets.csv"))
perm_13 = read.csv("perm_13.csv")
perm_17 = read.csv("perm_17.csv")

set.seed(1)
N = 200L
wd  = sample_wind_dir(N)$samples
ws  = sample_wind_speed(N)
ign = sample_ignition(landscape, N)
scenarios = data.frame(
  ign_row    = ign[, "row"],
  ign_col    = ign[, "col"],
  wind_speed = ws,
  wind_dir   = wd
)

cat("=== B=13, R=4 (N=200) ===\n")
t0 = proc.time()[3]
out = evaluate_all_plans(landscape  = landscape,
                         targets    = targets,
                         perm_13    = perm_13,
                         perm_17    = perm_17,
                         scenarios  = scenarios,
                         alpha      = 0.1,
                         R          = 4L)
cat(sprintf("Done in %.0fs\n\n", proc.time()[3] - t0))
saveRDS(out, "sensitivity_R4_N200.rds")
cat("Saved sensitivity_R4_N200.rds\n\n")

for (p in c("none", "perm_13", "perm_17", "reactive")) {
  d = out$damage[out$plan == p]
  d = d[!is.na(d)]
  rm = compute_risk_measures(d, alpha = 0.9, n_boot = 1000)
  cat(sprintf("  %-10s E=%.2f [%.2f,%.2f]  CVaR=%.2f [%.2f,%.2f]\n",
              p, rm$estimate[1], rm$ci_lower[1], rm$ci_upper[1],
                 rm$estimate[2], rm$ci_lower[2], rm$ci_upper[2]))
}

g = out$mip_gap[out$plan == "reactive"]
s = out$solver_status[out$plan == "reactive"]
cat(sprintf("\nReactive solver: %d/%d optimal, median gap=%.4f, max=%.4f\n",
            sum(s == 0, na.rm = TRUE), sum(!is.na(s)),
            median(g, na.rm = TRUE), max(g, na.rm = TRUE)))
