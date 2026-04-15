# ============================================================
# solve_reactive.R — Single-scenario reactive firebreak driver
# ============================================================
# Calls build_step4_dat() to write a wind-specific .dat, then
# runs step4_reactive.mod via AMPL + Gurobi and parses the
# reactive cells and final MIP gap from the output.
# ============================================================

source("build_step4_dat.R")


#' Solve one reactive-firebreak scenario
#'
#' Generates a wind-dependent .dat, launches AMPL with Gurobi
#' (timelim=120, mipgap=0.01), and parses the solution.
#'
#' @param landscape  Numeric matrix (21x21).
#' @param targets    Data frame with columns row, col, weight.
#' @param perm_cells Integer matrix (K x 2) — permanent cells from Step 3.
#' @param wind_speed Scalar >= 0. Revealed wind speed (km/h).
#' @param wind_dir   Scalar in [0, 2*pi). Direction wind blows FROM.
#' @param alpha      Scalar in (0, 1]. Compression exponent.
#' @param R          Integer >= 0. Reactive budget (default 4).
#' @param mod_file   Path to step4_reactive.mod.
#' @param timelim    Gurobi time limit in seconds (default 120).
#' @param mipgap     Gurobi relative MIP gap tolerance (default 0.01).
#' @param threads    Gurobi thread count per solve (default 1). Keep at 1
#'   when running under \code{mclapply} to avoid cores fighting for CPU.
#' @param ampl_bin   Path to AMPL binary. Defaults to \code{AMPL_BIN}
#'   environment variable.
#' @param gurobi_bin Path to Gurobi binary. Defaults to \code{GUROBI_BIN}
#'   environment variable.
#'
#' @return Named list with components:
#'   \describe{
#'     \item{damage}{Numeric objective value (weighted target reachability).}
#'     \item{status}{Integer solve result code (0 = optimal).}
#'     \item{mip_gap}{Numeric final MIP gap (0--1 scale), or NA.}
#'     \item{react_cells}{Character vector of reactively cleared cell names
#'       ("r_c" format), length 0..R.}
#'     \item{n_react}{Integer count of reactive cells.}
#'     \item{damaged}{Data frame with columns cell, weight, x for targets
#'       with x > 0.01, or NULL if none.}
#'   }
solve_reactive = function(landscape, targets, perm_cells,
                          wind_speed, wind_dir, alpha,
                          R = 4L,
                          mod_file = "step4_reactive.mod",
                          timelim = 120, mipgap = 0.01,
                          threads = 1L,
                          ampl_bin   = Sys.getenv("AMPL_BIN"),
                          gurobi_bin = Sys.getenv("GUROBI_BIN")) {

  if (ampl_bin   == "") stop("AMPL_BIN environment variable not set.")
  if (gurobi_bin == "") stop("GUROBI_BIN environment variable not set.")

  # --- 1. Write .dat -----------------------------------------
  dat_file = tempfile(fileext = ".dat")
  on.exit(unlink(dat_file), add = TRUE)

  build_step4_dat(landscape  = landscape,
                  targets    = targets,
                  perm_cells = perm_cells,
                  wind_speed = wind_speed,
                  wind_dir   = wind_dir,
                  alpha      = alpha,
                  R          = R,
                  outfile    = dat_file)

  # --- 2. Write temp .run ------------------------------------
  tmp_run = tempfile(fileext = ".run")
  on.exit(unlink(tmp_run), add = TRUE)

  run_lines = c(
    "reset;",
    sprintf("model %s;", mod_file),
    sprintf("data %s;",  dat_file),
    "",
    sprintf("option solver '%s';", gurobi_bin),
    sprintf("option gurobi_options 'timelim=%d mipgap=%g threads=%d';",
            timelim, mipgap, threads),
    "",
    "solve;",
    "",
    "# --- parseable output ------------------------------------",
    'printf "SUMMARY,%f,%d\\n", Total_Damage, solve_result_num;',
    "",
    "for {v in CELLS: react[v] > 0.5} {",
    '    printf "REACTED,%s\\n", v;',
    "}",
    "",
    "for {t in TARGETS: x[t] > 0.01} {",
    '    printf "DAMAGED,%s,%d,%f\\n", t, weight[t], x[t];',
    "}"
  )
  writeLines(run_lines, tmp_run)

  # --- 3. Run AMPL -------------------------------------------
  output = system2(ampl_bin, tmp_run, stdout = TRUE, stderr = TRUE)

  # --- 4. Parse SUMMARY --------------------------------------
  summary_line = grep("^SUMMARY,", output, value = TRUE)
  if (length(summary_line) == 0) {
    warning("No SUMMARY line — AMPL may have failed.")
    return(list(damage     = NA_real_,
                status     = NA_integer_,
                mip_gap    = NA_real_,
                react_cells = character(0),
                n_react    = 0L,
                damaged    = NULL))
  }

  parts  = strsplit(summary_line[1], ",")[[1]]
  damage = as.numeric(parts[2])
  status = as.integer(parts[3])

  # --- 5. Parse REACTED cells --------------------------------
  react_lines = grep("^REACTED,", output, value = TRUE)
  react_cells = sub("^REACTED,", "", react_lines)

  # --- 6. Parse DAMAGED targets ------------------------------
  dam_lines = grep("^DAMAGED,", output, value = TRUE)
  if (length(dam_lines) > 0) {
    dam_parts = strsplit(sub("^DAMAGED,", "", dam_lines), ",")
    damaged = data.frame(
      cell   = vapply(dam_parts, `[`, character(1), 1),
      weight = as.integer(vapply(dam_parts, `[`, character(1), 2)),
      x      = as.numeric(vapply(dam_parts, `[`, character(1), 3)),
      stringsAsFactors = FALSE
    )
  } else {
    damaged = NULL
  }

  # --- 7. Parse MIP gap from Gurobi output --------------------
  #    Gurobi 13 prints: "absmipgap=X, relmipgap=Y"
  gap_line = grep("relmipgap=", output, value = TRUE)
  if (length(gap_line) > 0) {
    gap_match = regmatches(gap_line[length(gap_line)],
                           regexpr("relmipgap=[0-9.eE+-]+",
                                   gap_line[length(gap_line)]))
    mip_gap = as.numeric(sub("relmipgap=", "", gap_match))
  } else {
    mip_gap = if (status == 0) 0 else NA_real_
  }

  # --- 8. Return ---------------------------------------------
  list(damage      = damage,
       status      = status,
       mip_gap     = mip_gap,
       react_cells = react_cells,
       n_react     = length(react_cells),
       damaged     = damaged)
}
