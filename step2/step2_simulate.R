# ============================================================
# step2_simulate.R — Validate MILP firebreak plans with
#                     Monte Carlo fire simulation
# ============================================================
# Reads optimal plans from step2_results.csv and
# step2_results_A.csv, de-duplicates by cleared-cell set,
# runs simulate_fire() K=200 times per plan under southern-
# front ignition (all row-21 cells) with no wind, and
# compares against baseline and current corridor.
# ============================================================

source("../.Rprofile")
source("../../Project data and resources/fire-simulator.r")

set.seed(324)

# --- Data ----------------------------------------------------
data_dir  <- "../../Project data and resources"
landscape <- as.matrix(read.table(file.path(data_dir, "landscape.csv"),
                                  sep = ",", header = FALSE))
targets   <- read.csv(file.path(data_dir, "targets.csv"))
# K: set before source(), or pass via command line, or defaults to 200
# Console:  K <- 500; source("step2_simulate.R")
# Terminal: Rscript step2_simulate.R 500
if (!exists("K")) {
  args <- commandArgs(trailingOnly = TRUE)
  K    <- if (length(args) >= 1) as.integer(args[1]) else 200L
}
nr        <- nrow(landscape)
nc        <- ncol(landscape)

# Southern front ignition: all row-21 vegetated cells
ignition_rc <- which(landscape[nr, ] != 0)
ignition_rc <- cbind(row = rep(nr, length(ignition_rc)),
                     col = ignition_rc)

# =============================================================
# Helpers
# =============================================================

#' Parse semicolon-separated "r_c" strings into (row, col) matrix
#'
#' @param cell_str Single string of semicolon-separated node names,
#'   e.g. "10_5;10_6;10_7"
#'
#' @return Integer matrix with columns named `row` and `col`
parse_cells <- function(cell_str) {
  tokens <- strsplit(cell_str, ";")[[1]]
  parts  <- strsplit(tokens, "_")
  rc     <- do.call(rbind, lapply(parts, as.integer))
  colnames(rc) <- c("row", "col")
  rc
}

#' Run K fire simulations and return mean damage with 95% CI
#'
#' @param landscape  21x21 matrix (firebreaks already applied as 0)
#' @param ignition   Two-column matrix (row, col) of ignition cells
#' @param targets    Data frame with row, col, weight columns
#' @param K          Number of replications
#'
#' @return Named list: mean_damage, se, ci_lo, ci_hi
run_simulation <- function(landscape, ignition, targets, K) {
  damages <- numeric(K)
  for (k in seq_len(K)) {
    result     <- simulate_fire(landscape, ignition,
                                wind_speed = 0, wind_dir = 0)
    damages[k] <- compute_damage(result$burned, targets)
  }
  mu <- mean(damages)
  se <- sd(damages) / sqrt(K)
  list(mean_damage = mu,
       se          = se,
       ci_lo       = mu - 1.96 * se,
       ci_hi       = mu + 1.96 * se)
}

# =============================================================
# 1. Collect unique firebreak plans from MILP results
# =============================================================

#' Read MILP result CSVs and extract unique non-empty plans
#'
#' @param files Character vector of CSV paths
#'
#' @return Data frame with one row per unique cleared-cell set, keeping
#'   the first occurrence's model/B/param metadata and MILP damage
collect_plans <- function(files) {
  all_res <- do.call(rbind, lapply(files, read.csv,
                                   stringsAsFactors = FALSE))

  # Keep feasible solves (0=optimal, 402=time limit with solution)
  all_res <- all_res[!is.na(all_res$status) &
                     all_res$status %in% c(0, 402), ]
  all_res <- all_res[!is.na(all_res$cleared_cells) &
                     nchar(trimws(all_res$cleared_cells)) > 0, ]

  # Normalise cell sets: sort cells within each string for dedup
  all_res$cells_key <- vapply(all_res$cleared_cells, function(s) {
    paste(sort(strsplit(s, ";")[[1]]), collapse = ";")
  }, character(1))

  # De-duplicate
  all_res <- all_res[!duplicated(all_res$cells_key), ]

  all_res$plan_label <- sprintf("%s_B%d_%s%s",
                                all_res$model, all_res$B,
                                all_res$param_name,
                                gsub("\\.", "",
                                     as.character(all_res$param_value)))
  all_res
}

plans <- collect_plans(c("step2_results.csv", "step2_results_A.csv"))
cat("Found", nrow(plans), "unique firebreak plans to simulate\n")

# =============================================================
# 2. Build plan list: baseline + corridor + MILP plans
# =============================================================
plan_list <- list()

# (a) Baseline — no firebreaks
plan_list[[1]] <- list(
  label = "baseline", break_rc = NULL,
  milp_damage = NA, model = "none", B = 0,
  param_name = NA, param_value = NA
)

# (b) Current corridor — row 10, cols 3-19 (17 centre cells)
corr_cols <- 3:19
plan_list[[2]] <- list(
  label = "corridor_row10",
  break_rc = cbind(row = rep(10, length(corr_cols)), col = corr_cols),
  milp_damage = NA, model = "corridor", B = length(corr_cols),
  param_name = NA, param_value = NA
)

# (c) MILP plans
for (i in seq_len(nrow(plans))) {
  plan_list[[length(plan_list) + 1]] <- list(
    label       = plans$plan_label[i],
    break_rc    = parse_cells(plans$cleared_cells[i]),
    milp_damage = plans$damage[i],
    model       = plans$model[i],
    B           = plans$B[i],
    param_name  = plans$param_name[i],
    param_value = plans$param_value[i]
  )
}

cat("Total plans to simulate:", length(plan_list), "\n")

# =============================================================
# 3. Simulate each plan
# =============================================================
sim_results <- vector("list", length(plan_list))

for (idx in seq_along(plan_list)) {
  pl <- plan_list[[idx]]
  cat(sprintf("[%d/%d] %s\n", idx, length(plan_list), pl$label))

  # Apply firebreaks to landscape
  ls_mod <- landscape
  if (!is.null(pl$break_rc)) {
    ls_mod <- set_firebreaks(landscape, pl$break_rc)
  }

  sim <- run_simulation(ls_mod, ignition_rc, targets, K)

  sim_results[[idx]] <- data.frame(
    plan_label  = pl$label,
    model       = pl$model,
    B           = pl$B,
    param_name  = ifelse(is.na(pl$param_name), "", pl$param_name),
    param_value = ifelse(is.na(pl$param_value), NA, pl$param_value),
    milp_damage = ifelse(is.na(pl$milp_damage), NA, pl$milp_damage),
    sim_mean    = sim$mean_damage,
    sim_se      = sim$se,
    sim_ci_lo   = sim$ci_lo,
    sim_ci_hi   = sim$ci_hi,
    n_cleared   = if (is.null(pl$break_rc)) 0L else nrow(pl$break_rc),
    stringsAsFactors = FALSE
  )
}

sim_results <- do.call(rbind, sim_results)

# =============================================================
# 4. Save
# =============================================================
outfile <- sprintf("step2_sim_results_K%d.csv", K)
write.csv(sim_results, outfile, row.names = FALSE)
cat(sprintf("\nDone. Saved %d rows to %s\n", nrow(sim_results), outfile))

# --- Quick summary -------------------------------------------
cat("\n=== Summary ===\n")
print(sim_results[, c("plan_label", "n_cleared", "milp_damage",
                       "sim_mean", "sim_ci_lo", "sim_ci_hi")],
      row.names = FALSE)
