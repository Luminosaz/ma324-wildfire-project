# ============================================================
# step2_solve_A_refine.R — Finer τ sweep for Model A only
# ============================================================
# First sweep showed τ=0~0.1 all burn, τ=0.15+ nothing burns.
# Need finer grid in 0.10–0.15 to find the transition.
# ============================================================

source("../.Rprofile")
source("step2_build_dat.R")

ampl_bin <- Sys.getenv("AMPL_BIN")
if (ampl_bin == "") stop("AMPL_BIN not set.")

B_vals   <- c(5, 8, 10, 13, 16, 20)
tau_vals <- c(0, 0.05, 0.10, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15, 0.20)

#' Run a single AMPL Model A solve
#'
#' @param B   Integer budget
#' @param tau Numeric threshold
#' @return One-row data frame with model, B, param_name, param_value,
#'   damage, n_cleared, cleared_cells, status
run_one_A <- function(B, tau) {
  tmp <- tempfile(fileext = ".run")
  writeLines(c(
    "reset;",
    "model step2_modelA.mod;",
    "data step2.dat;",
    sprintf("let B := %d;", B),
    sprintf("let tau := %g;", tau),
    "include step2_modelA.run;"
  ), tmp)

  output <- system2(ampl_bin, tmp, stdout = TRUE, stderr = TRUE)
  unlink(tmp)

  summary_line <- grep("^SUMMARY,", output, value = TRUE)
  if (length(summary_line) == 0) {
    warning(sprintf("No SUMMARY: B=%d tau=%g", B, tau))
    return(data.frame(
      model = "A", B = B, param_name = "tau",
      param_value = tau, damage = NA, n_cleared = NA,
      cleared_cells = NA, status = NA_integer_, stringsAsFactors = FALSE
    ))
  }

  parts   <- strsplit(summary_line, ",")[[1]]
  damage  <- as.numeric(parts[2])
  status  <- as.integer(parts[3])
  cleared <- sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE))

  data.frame(
    model = "A", B = B, param_name = "tau", param_value = tau,
    damage = damage, n_cleared = length(cleared),
    cleared_cells = paste(cleared, collapse = ";"),
    status = status, stringsAsFactors = FALSE
  )
}

# --- Run sweep -----------------------------------------------
results <- data.frame()
total <- length(B_vals) * length(tau_vals)
i <- 0

for (B in B_vals) {
  for (tau in tau_vals) {
    i <- i + 1
    cat(sprintf("[%d/%d] Model A  B=%d  tau=%g\n", i, total, B, tau))
    results <- rbind(results, run_one_A(B, tau))
  }
}

# --- Marginal value ------------------------------------------
results$marginal_value <- NA_real_
for (tau in tau_vals) {
  sub <- results[results$param_value == tau, ]
  sub <- sub[order(sub$B), ]
  for (j in 1:(nrow(sub) - 1)) {
    delta_B <- sub$B[j + 1] - sub$B[j]
    idx <- which(results$param_value == tau & results$B == sub$B[j])
    results$marginal_value[idx] <- (sub$damage[j] - sub$damage[j + 1]) / delta_B
  }
}

write.csv(results, "step2_results_A.csv", row.names = FALSE)
cat("\nDone. Saved", nrow(results), "rows to step2_results_A.csv\n")
