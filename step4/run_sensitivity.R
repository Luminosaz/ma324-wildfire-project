# ============================================================
# run_sensitivity.R — Step 4 Phase E: sensitivity analysis
# ============================================================
# 1. Sweep R ∈ {0, 2, 6, 8} with perm_13 + perm_17 (N=200).
#    R=4 already evaluated in final_N1000.rds — skip it.
# 2. (B=12, R=5) variant using perm_12 from perm_12.csv.
#
# Saves:  sensitivity_R_sweep.rds, sensitivity_B12_R5.rds
# ============================================================

source("../.Rprofile")

data_dir = "../../Project data and resources"
source(file.path(data_dir, "fire-simulator.r"))
source("../step3/step3_samplers.R")
source("../step3/step3_monte_carlo.R")
source("evaluate_all_plans.R")

# --- data -----------------------------------------------------
landscape = as.matrix(read.csv(file.path(data_dir, "landscape.csv"),
                               header = FALSE))
targets   = read.csv(file.path(data_dir, "targets.csv"))
perm_13   = read.csv("perm_13.csv")
perm_17   = read.csv("perm_17.csv")

alpha = 0.1
N     = 200L

# --- sample scenarios (single global seed) --------------------
set.seed(1)
wind_dir   = sample_wind_dir(N)$samples
wind_speed = sample_wind_speed(N)
ignition   = sample_ignition(landscape, N)

scenarios = data.frame(
  ign_row    = ignition[, "row"],
  ign_col    = ignition[, "col"],
  wind_speed = wind_speed,
  wind_dir   = wind_dir
)

cat(sprintf("Sampled %d scenarios for sensitivity runs\n\n", N))

# =============================================================
# 1. R sweep: R ∈ {0, 2, 6, 8}
# =============================================================
R_vals  = c(0L, 2L, 6L, 8L)
R_results = list()

for (R_val in R_vals) {
  cat(sprintf("=== R = %d ===\n", R_val))
  res = evaluate_all_plans(landscape  = landscape,
                           targets    = targets,
                           perm_13    = perm_13,
                           perm_17    = perm_17,
                           scenarios  = scenarios,
                           alpha      = alpha,
                           R          = R_val)
  res$R_budget = R_val
  R_results[[as.character(R_val)]] = res
}

sweep_df = do.call(rbind, R_results)
rownames(sweep_df) = NULL
saveRDS(sweep_df, "sensitivity_R_sweep.rds")
cat("\nSaved sensitivity_R_sweep.rds\n\n")

# =============================================================
# 2. (B=12, R=5) variant
# =============================================================
cat("=== B=12, R=5 ===\n")
perm_12 = read.csv("perm_12.csv")

res_b12 = evaluate_all_plans(landscape  = landscape,
                             targets    = targets,
                             perm_13    = perm_12,   # slot takes the permanent plan
                             perm_17    = perm_17,
                             scenarios  = scenarios,
                             alpha      = alpha,
                             R          = 5L)
res_b12$R_budget = 5L
res_b12$B_perm   = 12L
saveRDS(res_b12, "sensitivity_B12_R5.rds")
cat("Saved sensitivity_B12_R5.rds\n\n")

# =============================================================
# 3. Summary table: E[damage] and CVaR90 per plan per R
# =============================================================
cat("========================================\n")
cat("  Sensitivity summary\n")
cat("========================================\n\n")

summarise_plan = function(dmg, label) {
  dmg = dmg[!is.na(dmg)]
  if (length(dmg) == 0) return(data.frame(plan = label,
                                           E_damage = NA, CVaR90 = NA))
  rm = compute_risk_measures(dmg, alpha = 0.9, n_boot = 2000)
  data.frame(plan     = label,
             E_damage = round(rm$estimate[1], 2),
             CVaR90   = round(rm$estimate[2], 2))
}

# --- R sweep summary ------------------------------------------
cat("--- R sweep (N=200, B=13) ---\n")
for (R_val in R_vals) {
  cat(sprintf("\nR = %d:\n", R_val))
  sub = sweep_df[sweep_df$R_budget == R_val, ]
  tbl = do.call(rbind, lapply(unique(sub$plan), function(p) {
    summarise_plan(sub$damage[sub$plan == p], p)
  }))
  print(tbl, row.names = FALSE)
}

# --- B=12, R=5 summary ----------------------------------------
cat("\n--- B=12, R=5 (N=200) ---\n")
tbl_b12 = do.call(rbind, lapply(unique(res_b12$plan), function(p) {
  summarise_plan(res_b12$damage[res_b12$plan == p], p)
}))
print(tbl_b12, row.names = FALSE)

cat("\nDone.\n")
