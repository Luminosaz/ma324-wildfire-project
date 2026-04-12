# ============================================================
# step2_solve.R — Sweep B × τ (Model A) and B × α (Model B)
# ============================================================

source("../.Rprofile")
source("step2_build_dat.R")

ampl_bin = Sys.getenv("AMPL_BIN")
if (ampl_bin == "") stop("AMPL_BIN environment variable not set.")

# --- Sweep grids ---------------------------------------------
B_vals    = c(5, 8, 10, 13, 16, 20)
tau_vals  = c(0, 0.05, 0.10, 0.15, 0.20)
alpha_vals = c(0.05, 0.1, 0.15, 0.2, 0.3)

# =============================================================
# Helper: run one AMPL solve and parse output
# =============================================================

#' Run a single AMPL firebreak solve
#'
#' @param model      Character "A" or "B"
#' @param B          Integer budget
#' @param param_name Character "tau" or "alpha"
#' @param param_val  Numeric value for the swept parameter
#'
#' @return A one-row data frame with columns: model, B, param_name,
#'   param_value, damage, n_cleared, cleared_cells, status
run_one = function(model, B, param_name, param_val) {
  mod_file = sprintf("step2_model%s.mod", model)
  run_file = sprintf("step2_model%s.run", model)

  tmp = tempfile(fileext = ".run")
  lines = c(
    "reset;",
    sprintf("model %s;", mod_file),
    "data step2.dat;",
    sprintf("let B := %d;", B),
    sprintf("let %s := %g;", param_name, param_val),
    sprintf("include %s;", run_file)
  )
  writeLines(lines, tmp)

  output = system2(ampl_bin, tmp, stdout = TRUE, stderr = TRUE)
  unlink(tmp)

  # Parse SUMMARY line
  summary_line = grep("^SUMMARY,", output, value = TRUE)
  if (length(summary_line) == 0) {
    warning(sprintf("No SUMMARY: model=%s B=%d %s=%g", model, B, param_name, param_val))
    return(data.frame(
      model = model, B = B, param_name = param_name,
      param_value = param_val, damage = NA, n_cleared = NA,
      cleared_cells = NA, status = NA_integer_, stringsAsFactors = FALSE
    ))
  }

  parts = strsplit(summary_line, ",")[[1]]
  damage = as.numeric(parts[2])
  status = as.integer(parts[3])    # solve_result_num: 0=optimal

  # Parse CLEARED lines
  cleared = sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE))

  data.frame(
    model        = model,
    B            = B,
    param_name   = param_name,
    param_value  = param_val,
    damage       = damage,
    n_cleared    = length(cleared),
    cleared_cells = paste(cleared, collapse = ";"),
    status       = status,
    stringsAsFactors = FALSE
  )
}

# =============================================================
# 2. Run all combos
# =============================================================
results = data.frame()
total  = length(B_vals) * (length(tau_vals) + length(alpha_vals))
i      = 0

# Model A: B × τ
for (B in B_vals) {
  for (tau in tau_vals) {
    i = i + 1
    cat(sprintf("[%d/%d] Model A  B=%d  tau=%g\n", i, total, B, tau))
    results = rbind(results, run_one("A", B, "tau", tau))
  }
}

# Model B: B × α
for (B in B_vals) {
  for (alpha in alpha_vals) {
    i = i + 1
    cat(sprintf("[%d/%d] Model B  B=%d  alpha=%g\n", i, total, B, alpha))
    results = rbind(results, run_one("B", B, "alpha", alpha))
  }
}

# =============================================================
# 3. Compute marginal value of one extra cleared cell
# =============================================================

#' Compute marginal damage reduction for consecutive budget levels
#'
#' @param df  Data frame with at least columns: model, param_value, B, damage
#'
#' @return The input data frame with an additional column `marginal_value`,
#'   defined as damage(B_current) - damage(B_next) for consecutive B values
#'   within each (model, param_value) group. NA for the largest B in each group.
add_marginal_value = function(df) {
  df$marginal_value = NA_real_
  groups = split(df, list(df$model, df$param_value))

  for (g in groups) {
    if (nrow(g) < 2) next
    g = g[order(g$B), ]
    for (j in 1:(nrow(g) - 1)) {
      delta_B = g$B[j + 1] - g$B[j]
      idx = which(df$model == g$model[j] &
                   df$param_value == g$param_value[j] &
                   df$B == g$B[j])
      df$marginal_value[idx] = (g$damage[j] - g$damage[j + 1]) / delta_B
    }
  }
  df
}

results = add_marginal_value(results)

# =============================================================
# 4. Save
# =============================================================
write.csv(results, "step2_results.csv", row.names = FALSE)
cat("\nDone. Saved", nrow(results), "rows to step2_results.csv\n")
