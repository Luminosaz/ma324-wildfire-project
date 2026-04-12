# ============================================================
# test_modelB.R — Sanity check: run Model B once (B=13, alpha=0.5)
# ============================================================

# --- 1. Write temp .run wrapper ------------------------------
tmp_run <- tempfile(fileext = ".run")
writeLines(c(
  "reset;",
  "model step2_modelB.mod;",
  "data step2.dat;",
  "let B := 13;",
  "let alpha := 0.5;",
  "include step2_modelB.run;"
), tmp_run)

# --- 2. Run AMPL --------------------------------------------
output <- system2(Sys.getenv("AMPL_BIN"), tmp_run, stdout = TRUE, stderr = TRUE)
cat("=== Raw AMPL output ===\n")
cat(output, sep = "\n")

# --- 3. Parse results ----------------------------------------
summary_line  <- grep("^SUMMARY,", output, value = TRUE)
cleared_lines <- grep("^CLEARED,", output, value = TRUE)
damaged_lines <- grep("^DAMAGED,", output, value = TRUE)

if (length(summary_line) == 0) {
  stop("No SUMMARY line found — AMPL may have failed. Check output above.")
}

parts <- strsplit(summary_line, ",")[[1]]
cat("\n=== Results ===\n")
cat("Objective (total damage):", parts[2], "\n")
cat("Solver message:", parts[3], "\n")
cat("Cells cleared:", length(cleared_lines), "\n")
cat("  ", sub("^CLEARED,", "", cleared_lines), "\n")
cat("Targets damaged:", length(damaged_lines), "\n")
for (d in damaged_lines) {
  dp <- strsplit(sub("^DAMAGED,", "", d), ",")[[1]]
  cat("  ", dp[1], "(weight=", dp[2], ", reachability=", dp[3], ")\n")
}
