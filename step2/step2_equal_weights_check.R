# =============================================================================
# Equal-weights sensitivity (brief §4.2.4 Q4): Model B, B=13, alpha=0.1, w=5:5
# =============================================================================
# Re-solves Model B with settlement and wetland weights both set to 5, to
# test whether the L-shape geometry is weight-driven or geometry-driven.
# Runs 120s (sweep-matching) and 600s (extended) solves.

source("../.Rprofile")
source("step2_build_dat.R")

targets_eq = targets
targets_eq$weight = 5
generate_step2_dat(edges, targets_eq, outfile = "step2_equal.dat",
                   landscape = landscape)

ampl_bin   = Sys.getenv("AMPL_BIN")
gurobi_bin = Sys.getenv("GUROBI_BIN")

#' Solve Model B at B=13, alpha=0.1 on a given .dat with a given time limit.
solve_B_13 = function(dat_file, timelim) {
  tmp = tempfile(fileext = ".run")
  writeLines(c(
    "reset;",
    "model step2_modelB.mod;",
    sprintf("data %s;", dat_file),
    "let B := 13;",
    "let alpha := 0.1;",
    sprintf("option solver '%s';", gurobi_bin),
    sprintf("option gurobi_options 'mipgap=0.01 timelim=%d';", timelim),
    "solve;",
    'printf "SUMMARY,%f,%d\\n", Total_Damage, solve_result_num;',
    'for {v in CELLS: clear[v] > 0.5} { printf "CLEARED,%s\\n", v; }',
    'for {t in TARGETS: x[t] > 0.01} { printf "DAMAGED,%s,%d,%f\\n", t, weight[t], x[t]; }'
  ), tmp)
  output = system2(ampl_bin, tmp, stdout = TRUE, stderr = TRUE)
  unlink(tmp)
  list(
    summary = grep("^SUMMARY,", output, value = TRUE),
    cleared = sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE)),
    damaged = grep("^DAMAGED,", output, value = TRUE),
    gap     = grep("relmipgap=", output, value = TRUE),
    full    = output
  )
}

cat("=== Solve 1: 120s (sweep-matching) ===\n")
r1 = solve_B_13("step2_equal.dat", 120)
cat("SUMMARY:", r1$summary, "\n")
cat("n_cleared:", length(r1$cleared), "\n")
cat("Cleared:", paste(r1$cleared, collapse = ", "), "\n")
cat("Gap:", r1$gap, "\n\n")

cat("=== Solve 2: 600s (extended) ===\n")
r2 = solve_B_13("step2_equal.dat", 600)
cat("SUMMARY:", r2$summary, "\n")
cat("n_cleared:", length(r2$cleared), "\n")
cat("Cleared:", paste(r2$cleared, collapse = ", "), "\n")
cat("Gap:", r2$gap, "\n\n")

# Persist full results
writeLines(c(
  "# Equal-weights sensitivity check (brief §4.2.4 Q4)",
  "# Model B, B=13, alpha=0.1, settlement weight = wetland weight = 5",
  sprintf("# Generated: %s", Sys.time()),
  "",
  "## Original weights (from targets.csv): settlement=10, wetland=5",
  "## Equal weights: settlement=5, wetland=5",
  "",
  "## Solve 1: 120s time limit (matches main sweep in step2_solve.R)",
  r1$summary, r1$gap,
  paste("CLEARED:", paste(r1$cleared, collapse = ", ")),
  r1$damaged,
  "",
  "## Solve 2: 600s time limit (extended check for incumbent stability)",
  r2$summary, r2$gap,
  paste("CLEARED:", paste(r2$cleared, collapse = ", ")),
  r2$damaged,
  "",
  "## Interpretation",
  "# L-shape reference (original 10:5 weights at B=13, alpha=0.1):",
  "#   col 15 rows 1-7 + row 7 cols 15-21 (settlement-focused)",
  sprintf("# L-shape match at 120s: %s",
          if (setequal(r1$cleared,
                       c(paste0(1:7, "_15"), paste0("7_", 16:21))))
            "YES (same geometry)" else "NO (different cells)"),
  sprintf("# L-shape match at 600s: %s",
          if (setequal(r2$cleared,
                       c(paste0(1:7, "_15"), paste0("7_", 16:21))))
            "YES (same geometry)" else "NO (different cells)")
), "step2_equal_weights_results.txt")
cat("Saved to step2/step2_equal_weights_results.txt\n")
