# ============================================================
# evaluate_all_plans.R — Step 4 Phase D: batch evaluation
# ============================================================
# Evaluates four firebreak plans under the same pre-sampled
# scenarios.  The caller owns set.seed(1) for scenario sampling.
# Within each scenario, set.seed(1000 + i) aligns the RNG
# starting state before each simulate_fire.  This does NOT create
# exact pathwise coupling — different firebreak layouts change
# the order and number of runif calls inside simulate_fire, so
# the four streams drift apart after a few spread decisions —
# but it reduces simulator noise when comparing plans at the
# scenario level and makes the batch fully reproducible.
#
# Plans:
#   1. none     — no firebreaks
#   2. perm_13  — B=13 permanent cells only
#   3. perm_17  — B+R=17 permanent cells (all budget pre-season)
#   4. reactive — B=13 permanent + up to R reactive (MILP per scenario)
# ============================================================

source("solve_reactive.R")     # also pulls in build_step4_dat.R, edges_wind.R


## ---------------------------------------------------------------
## Helper: parse "r_c" node strings to a two-column integer matrix
## ---------------------------------------------------------------
parse_nodes = function(nodes) {
  if (length(nodes) == 0) return(matrix(integer(0), ncol = 2,
                                        dimnames = list(NULL, c("row", "col"))))
  parts = strsplit(nodes, "_")
  rc    = do.call(rbind, lapply(parts, as.integer))
  colnames(rc) = c("row", "col")
  rc
}


#' Evaluate four firebreak plans under shared scenarios (CRN)
#'
#' For each row of \code{scenarios}, simulates fire under four plans
#' with the same ignition, wind, and internal RNG state.  Plan 4
#' (reactive) solves a per-scenario MILP; failures are caught by
#' \code{tryCatch} so a single bad solve does not kill the batch.
#'
#' @param landscape  Numeric matrix (21x21).
#' @param targets    Data frame with columns row, col, weight.
#' @param perm_13    Permanent cells (B=13). Matrix or data frame with
#'   row, col columns.
#' @param perm_17    Permanent cells (B+R=17). Matrix or data frame with
#'   row, col columns.
#' @param scenarios  Data frame with columns \code{ign_row}, \code{ign_col},
#'   \code{wind_speed}, \code{wind_dir}.  Pre-sampled by the caller under
#'   a single \code{set.seed}.
#' @param alpha      Scalar in (0, 1]. Compression exponent for reactive MILP.
#' @param R          Integer reactive budget (default 4).
#' @param mod_file   Path to step4_reactive.mod (default
#'   \code{"step4_reactive.mod"}).
#' @param timelim    Gurobi time limit in seconds (default 120).
#' @param mipgap     Gurobi relative MIP gap tolerance (default 0.01).
#' @param ampl_bin   Path to AMPL binary (default \code{AMPL_BIN} env var).
#' @param gurobi_bin Path to Gurobi binary (default \code{GUROBI_BIN} env var).
#'
#' @return A data frame (long format, N x 4 rows) with columns:
#'   \code{scenario_id}, \code{plan}, \code{damage},
#'   \code{react_cells} (list-column, NA for plans 1--3),
#'   \code{solve_time}, \code{mip_gap}, \code{solver_status}.
evaluate_all_plans = function(landscape, targets, perm_13, perm_17,
                              scenarios, alpha, R = 4L,
                              mod_file   = "step4_reactive.mod",
                              timelim    = 120,
                              mipgap     = 0.01,
                              ampl_bin   = Sys.getenv("AMPL_BIN"),
                              gurobi_bin = Sys.getenv("GUROBI_BIN")) {

  N = nrow(scenarios)
  plan_names = c("none", "perm_13", "perm_17", "reactive")
  n_plans    = length(plan_names)
  total_rows = N * n_plans

  # --- normalise perm_cells to matrix -------------------------
  if (is.data.frame(perm_13)) perm_13 = as.matrix(perm_13[, c("row", "col")])
  if (is.data.frame(perm_17)) perm_17 = as.matrix(perm_17[, c("row", "col")])

  # --- pre-build landscapes for plans 1--3 --------------------
  land_none = landscape
  land_13   = set_firebreaks(landscape, perm_13)
  land_17   = set_firebreaks(landscape, perm_17)

  # --- pre-allocate output vectors ----------------------------
  out_scen   = integer(total_rows)
  out_plan   = character(total_rows)
  out_damage = numeric(total_rows)
  out_react  = vector("list", total_rows)
  out_stime  = rep(NA_real_,    total_rows)
  out_gap    = rep(NA_real_,    total_rows)
  out_status = rep(NA_integer_, total_rows)

  cat(sprintf("Evaluating %d scenarios x %d plans ...\n", N, n_plans))
  t_start = proc.time()[3]

  for (i in seq_len(N)) {

    ir  = scenarios$ign_row[i]
    ic  = scenarios$ign_col[i]
    ign = matrix(c(ir, ic), nrow = 1,
                 dimnames = list(NULL, c("row", "col")))
    V   = scenarios$wind_speed[i]
    dir = scenarios$wind_dir[i]
    sim_seed = 1000L + i

    # --- Plan 1: none ----------------------------------------
    idx = (i - 1) * n_plans + 1
    out_scen[idx]   = i
    out_plan[idx]   = "none"
    out_react[[idx]] = NA_character_
    if (land_none[ir, ic] == 0) {
      out_damage[idx] = 0
    } else {
      set.seed(sim_seed)
      res1 = simulate_fire(land_none, ign, V, dir)
      out_damage[idx] = compute_damage(res1$burned, targets)
    }

    # --- Plan 2: perm_13 ------------------------------------
    idx = (i - 1) * n_plans + 2
    out_scen[idx]   = i
    out_plan[idx]   = "perm_13"
    out_react[[idx]] = NA_character_
    if (land_13[ir, ic] == 0) {
      out_damage[idx] = 0
    } else {
      set.seed(sim_seed)
      res2 = simulate_fire(land_13, ign, V, dir)
      out_damage[idx] = compute_damage(res2$burned, targets)
    }

    # --- Plan 3: perm_17 ------------------------------------
    idx = (i - 1) * n_plans + 3
    out_scen[idx]   = i
    out_plan[idx]   = "perm_17"
    out_react[[idx]] = NA_character_
    if (land_17[ir, ic] == 0) {
      out_damage[idx] = 0
    } else {
      set.seed(sim_seed)
      res3 = simulate_fire(land_17, ign, V, dir)
      out_damage[idx] = compute_damage(res3$burned, targets)
    }

    # --- Plan 4: reactive ------------------------------------
    idx = (i - 1) * n_plans + 4
    out_scen[idx] = i
    out_plan[idx] = "reactive"

    if (land_13[ir, ic] == 0) {
      # ignition lands on a permanent firebreak — no fire, skip MILP
      out_damage[idx]  = 0
      out_stime[idx]   = 0
      out_gap[idx]     = NA_real_
      out_status[idx]  = NA_integer_
      out_react[[idx]] = NA_character_
    } else {
      t0 = proc.time()[3]

      react_result = tryCatch(
        solve_reactive(landscape   = landscape,
                       targets     = targets,
                       perm_cells  = perm_13,
                       wind_speed  = V,
                       wind_dir    = dir,
                       alpha       = alpha,
                       R           = R,
                       mod_file    = mod_file,
                       timelim     = timelim,
                       mipgap      = mipgap,
                       ampl_bin    = ampl_bin,
                       gurobi_bin  = gurobi_bin),
        error = function(e) {
          warning(sprintf("Scenario %d: solve_reactive failed — %s", i, e$message))
          NULL
        }
      )

      out_stime[idx] = proc.time()[3] - t0

      if (!is.null(react_result) && !is.na(react_result$damage)) {
        out_gap[idx]    = react_result$mip_gap
        out_status[idx] = react_result$status
        out_react[[idx]] = react_result$react_cells

        # build landscape with permanent + reactive cells
        react_rc   = parse_nodes(react_result$react_cells)
        land_react = set_firebreaks(land_13, react_rc)

        if (land_react[ir, ic] == 0) {
          out_damage[idx] = 0
        } else {
          set.seed(sim_seed)
          res4 = simulate_fire(land_react, ign, V, dir)
          out_damage[idx] = compute_damage(res4$burned, targets)
        }
      } else {
        out_gap[idx]     = NA_real_
        out_status[idx]  = NA_integer_
        out_react[[idx]] = NA_character_
        out_damage[idx]  = NA_real_
      }
    }

    # --- progress --------------------------------------------
    if (i %% 50 == 0 || i == N) {
      cat(sprintf("  scenario %d / %d  (elapsed %.0fs)\n",
                  i, N, proc.time()[3] - t_start))
    }
  }

  # --- assemble data frame ------------------------------------
  out = data.frame(
    scenario_id   = out_scen,
    plan          = out_plan,
    damage        = out_damage,
    solve_time    = out_stime,
    mip_gap       = out_gap,
    solver_status = out_status,
    stringsAsFactors = FALSE
  )
  out$react_cells = out_react

  cat(sprintf("Done. %d rows (%d scenarios x %d plans)\n",
              nrow(out), N, n_plans))

  out
}
