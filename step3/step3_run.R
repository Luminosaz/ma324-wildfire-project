## ============================================================
## step3_run.R
## Runner script for Step 3: Monte Carlo evaluation of
## firebreak plans under uncertain wind and ignition.
## ============================================================

set.seed(1)

## ----------------------------------------------------------
## 0. Source dependencies
## ----------------------------------------------------------
source("step3_samplers.R")
source("step3_monte_carlo.R")

data_dir = "../../Project data and resources"
source(file.path(data_dir, "fire-simulator.r"))

## ----------------------------------------------------------
## 1. Load data
## ----------------------------------------------------------
landscape = as.matrix(read.csv(file.path(data_dir, "landscape.csv"),
                               header = FALSE))
targets   = read.csv(file.path(data_dir, "targets.csv"))

cat("Landscape:", nrow(landscape), "x", ncol(landscape), "\n")
cat("Targets:  ", nrow(targets), "rows\n")
cat("Vegetated:", sum(landscape != 0), "cells\n\n")

## ----------------------------------------------------------
## 2. Parse the proposed plan from Step 2 results
## ----------------------------------------------------------
step2 = read.csv("../step2/step2_results.csv")

## filter for model B, B=13, alpha=0.1
row_idx = which(step2$model      == "B"     &
                step2$B          == 13      &
                step2$param_name == "alpha" &
                step2$param_value == 0.1)

if (length(row_idx) != 1) {
  stop("Expected exactly 1 matching row in step2_results.csv, found ", length(row_idx))
}

## parse "10_5;10_6;10_7" into a 2-column matrix
parse_cells = function(s) {
  tokens = strsplit(s, ";")[[1]]
  parts  = strsplit(tokens, "_")
  rc     = do.call(rbind, lapply(parts, function(p) as.integer(p)))
  colnames(rc) = c("row", "col")
  rc
}

proposed_breaks = parse_cells(step2$cleared_cells[row_idx])
cat("Proposed plan (model B, B=13, alpha=0.1):\n")
cat("  Cells cleared:", nrow(proposed_breaks), "\n")
cat("  Locations:    ", step2$cleared_cells[row_idx], "\n\n")

## ----------------------------------------------------------
## 3. Define the three plans
## ----------------------------------------------------------
## Corridor: row 10, cols 3-19 (17 cells, current solution)
corridor_breaks = cbind(row = rep(10L, 17), col = 3L:19L)

plans = list(
  baseline = NULL,
  corridor = corridor_breaks,
  proposed = proposed_breaks
)

cat("Plans defined:\n")
for (nm in names(plans)) {
  n_cells = if (is.null(plans[[nm]])) 0 else nrow(plans[[nm]])
  cat(sprintf("  %-10s %2d cells\n", nm, n_cells))
}
cat("\n")

## ----------------------------------------------------------
## 4. Run Monte Carlo (default N = 1000, override before source)
## ----------------------------------------------------------
if (!exists("N")) N = 1000
cat("Running Monte Carlo with N =", N, "scenarios ...\n\n")

t0      = proc.time()
results = run_mc_evaluation(landscape, targets, N, plans)
elapsed = (proc.time() - t0)["elapsed"]

cat(sprintf("\nTotal elapsed time: %.1f seconds\n\n", elapsed))

## ----------------------------------------------------------
## 5. Compute risk measures for each plan
## ----------------------------------------------------------
risk_list = lapply(names(results), function(nm) {
  rm          = compute_risk_measures(results[[nm]]$damage)
  rm$plan     = nm
  rm
})

risk_table = do.call(rbind, risk_list)
risk_table = risk_table[, c("plan", "measure", "estimate", "ci_lower", "ci_upper")]

## ----------------------------------------------------------
## 6. Print comparison table
## ----------------------------------------------------------
cat("============================================================\n")
cat("  RISK MEASURE COMPARISON\n")
cat("============================================================\n\n")

for (m in unique(risk_table$measure)) {
  sub = risk_table[risk_table$measure == m, ]
  cat(m, ":\n")
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("  %-10s  %7.2f   [%7.2f, %7.2f]\n",
                sub$plan[j], sub$estimate[j], sub$ci_lower[j], sub$ci_upper[j]))
  }
  cat("\n")
}

## ----------------------------------------------------------
## 7. Save results
## ----------------------------------------------------------
## per-scenario results (long format)
all_scenarios = do.call(rbind, lapply(names(results), function(nm) {
  df      = results[[nm]]
  df$plan = nm
  df
}))
write.csv(all_scenarios, sprintf("step3_scenarios_N%d.csv", N), row.names = FALSE)

## risk summary
write.csv(risk_table, sprintf("step3_risk_measures_N%d.csv", N), row.names = FALSE)

cat("Saved:\n")
cat(sprintf("  step3_scenarios_N%d.csv      — per-scenario damage for all plans\n", N))
cat(sprintf("  step3_risk_measures_N%d.csv  — mean & CVaR with 95%% CIs\n", N))
