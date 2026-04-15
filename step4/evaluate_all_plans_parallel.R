# ============================================================
# evaluate_all_plans_parallel.R — Step 4 Phase D batch evaluation
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
# Scenarios are dispatched in parallel via parallel::mclapply
# (fork-based; works on macOS / Linux).  Per-scenario seeding
# gives identical results regardless of core count.
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
#' with the same ignition, wind, and per-scenario RNG seed.  Plan 4
#' (reactive) solves a per-scenario MILP; failures are caught by
#' \code{tryCatch} so a single bad solve does not kill the batch.
#' Scenarios are dispatched in parallel with \code{parallel::mclapply}.
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
#' @param threads    Gurobi thread count per solve (default 1). Passed to
#'   \code{solve_reactive}; keep at 1 under parallel dispatch.
#' @param n_cores    Integer number of cores for \code{mclapply}.
#'   Default \code{parallel::detectCores() - 1L} (leave one for the OS).
#'   Set to 1 for sequential execution (useful for debugging).
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
                              threads    = 1L,
                              ampl_bin   = Sys.getenv("AMPL_BIN"),
                              gurobi_bin = Sys.getenv("GUROBI_BIN"),
                              n_cores    = max(1L, parallel::detectCores() - 1L)) {

  N = nrow(scenarios)
  plan_names = c("none", "perm_13", "perm_17", "reactive")

  # --- normalise perm_cells to matrix -------------------------
  if (is.data.frame(perm_13)) perm_13 = as.matrix(perm_13[, c("row", "col")])
  if (is.data.frame(perm_17)) perm_17 = as.matrix(perm_17[, c("row", "col")])

  # --- pre-build landscapes for plans 1--3 --------------------
  land_none = landscape
  land_13   = set_firebreaks(landscape, perm_13)
  land_17   = set_firebreaks(landscape, perm_17)

  # --- worker: evaluate one scenario --------------------------
  run_one = function(i) {

    ir  = scenarios$ign_row[i]
    ic  = scenarios$ign_col[i]
    ign = matrix(c(ir, ic), nrow = 1,
                 dimnames = list(NULL, c("row", "col")))
    V   = scenarios$wind_speed[i]
    dir = scenarios$wind_dir[i]
    sim_seed = 1000L + i

    # per-scenario output
    dmg    = numeric(4)
    rcells = vector("list", 4)
    stime  = rep(NA_real_, 4)
    gap    = rep(NA_real_, 4)
    status = rep(NA_integer_, 4)

    # Plan 1: none
    rcells[[1]] = NA_character_
    if (land_none[ir, ic] == 0) {
      dmg[1] = 0
    } else {
      set.seed(sim_seed)
      res = simulate_fire(land_none, ign, V, dir)
      dmg[1] = compute_damage(res$burned, targets)
    }

    # Plan 2: perm_13
    rcells[[2]] = NA_character_
    if (land_13[ir, ic] == 0) {
      dmg[2] = 0
    } else {
      set.seed(sim_seed)
      res = simulate_fire(land_13, ign, V, dir)
      dmg[2] = compute_damage(res$burned, targets)
    }

    # Plan 3: perm_17
    rcells[[3]] = NA_character_
    if (land_17[ir, ic] == 0) {
      dmg[3] = 0
    } else {
      set.seed(sim_seed)
      res = simulate_fire(land_17, ign, V, dir)
      dmg[3] = compute_damage(res$burned, targets)
    }

    # Plan 4: reactive
    if (land_13[ir, ic] == 0) {
      dmg[4]      = 0
      stime[4]    = 0
      rcells[[4]] = NA_character_
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
                       gurobi_bin  = gurobi_bin,
                       threads     = threads),
        error = function(e) {
          warning(sprintf("Scenario %d: solve_reactive failed — %s",
                          i, e$message))
          NULL
        }
      )

      stime[4] = proc.time()[3] - t0

      if (!is.null(react_result) && !is.na(react_result$damage)) {
        gap[4]      = react_result$mip_gap
        status[4]   = react_result$status
        rcells[[4]] = react_result$react_cells

        react_rc   = parse_nodes(react_result$react_cells)
        land_react = set_firebreaks(land_13, react_rc)

        if (land_react[ir, ic] == 0) {
          dmg[4] = 0
        } else {
          set.seed(sim_seed)
          res = simulate_fire(land_react, ign, V, dir)
          dmg[4] = compute_damage(res$burned, targets)
        }
      } else {
        gap[4]      = NA_real_
        status[4]   = NA_integer_
        rcells[[4]] = NA_character_
        dmg[4]      = NA_real_
      }
    }

    list(scenario_id   = rep(i, 4),
         plan          = plan_names,
         damage        = dmg,
         react_cells   = rcells,
         solve_time    = stime,
         mip_gap       = gap,
         solver_status = status)
  }

  # --- dispatch -----------------------------------------------
  cat(sprintf("Evaluating %d scenarios x 4 plans on %d core%s ...\n",
              N, n_cores, if (n_cores > 1) "s" else ""))
  t_start = proc.time()[3]

  results = parallel::mclapply(seq_len(N), run_one,
                               mc.cores       = n_cores,
                               mc.preschedule = FALSE)

  elapsed = proc.time()[3] - t_start
  cat(sprintf("Done. %.0fs elapsed (%.1fs / scenario avg)\n",
              elapsed, elapsed / N))

  # --- check for mclapply errors ------------------------------
  bad = vapply(results, inherits, logical(1), "try-error")
  if (any(bad)) {
    warning(sprintf("%d scenario(s) returned try-error; replaced with NA.",
                    sum(bad)))
    na_row = list(scenario_id = rep(NA_integer_, 4),
                  plan        = plan_names,
                  damage      = rep(NA_real_, 4),
                  react_cells = as.list(rep(NA_character_, 4)),
                  solve_time  = rep(NA_real_, 4),
                  mip_gap     = rep(NA_real_, 4),
                  solver_status = rep(NA_integer_, 4))
    for (j in which(bad)) {
      na_row$scenario_id = rep(j, 4)
      results[[j]] = na_row
    }
  }

  # --- combine ------------------------------------------------
  out = data.frame(
    scenario_id   = unlist(lapply(results, `[[`, "scenario_id")),
    plan          = unlist(lapply(results, `[[`, "plan")),
    damage        = unlist(lapply(results, `[[`, "damage")),
    solve_time    = unlist(lapply(results, `[[`, "solve_time")),
    mip_gap       = unlist(lapply(results, `[[`, "mip_gap")),
    solver_status = unlist(lapply(results, `[[`, "solver_status")),
    stringsAsFactors = FALSE
  )
  out$react_cells = do.call(c, lapply(results, `[[`, "react_cells"))

  cat(sprintf("%d rows (%d scenarios x %d plans)\n",
              nrow(out), N, length(plan_names)))

  out
}
