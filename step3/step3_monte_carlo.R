# ============================================================
# step3_monte_carlo.R — Monte Carlo evaluation of firebreak
#   plans under uncertain wind and ignition.
#
# Dependencies:
#   source("step3_samplers.R")
#   source("../../Project data and resources/fire-simulator.r")
# ============================================================


#' Monte Carlo evaluation of firebreak plans
#'
#' Pre-samples N scenarios (wind speed, wind direction, ignition point),
#' then evaluates every plan against the same random draws (common
#' random numbers) so that plan comparisons are fair.
#'
#' @param landscape Integer matrix (21x21) from landscape.csv.
#' @param targets Data frame with columns row, col, weight.
#' @param N Integer, number of Monte Carlo scenarios.
#' @param plans Named list of firebreak specifications. Each element is
#'   either NULL (no firebreaks) or an n x 2 integer matrix with
#'   columns row, col giving the cells to clear.
#' @return Named list (one entry per plan) of data frames with columns
#'   scenario, damage, wind_speed, wind_dir, ign_row, ign_col.
run_mc_evaluation = function(landscape, targets, N, plans) {

  ## 1. Pre-sample all N scenarios (shared across plans = CRN)
  cat("Sampling", N, "scenarios ...\n")
  wind_dir   = sample_wind_dir(N)$samples
  wind_speed = sample_wind_speed(N)
  ignition   = sample_ignition(landscape, N)   # N x 2 matrix

  ## 2. Evaluate each plan
  results = list()

  for (plan_name in names(plans)) {
    cat("\n--- Plan:", plan_name, "---\n")

    ls_plan = landscape
    breaks  = plans[[plan_name]]
    if (!is.null(breaks)) {
      ls_plan = set_firebreaks(ls_plan, breaks)
    }

    damages = numeric(N)

    for (i in seq_len(N)) {
      ir = ignition[i, "row"]
      ic = ignition[i, "col"]

      if (ls_plan[ir, ic] == 0) {
        ## ignition landed on a cell cleared by this plan — no fire
        damages[i] = 0
      } else {
        result     = simulate_fire(ls_plan, ignition[i, , drop = FALSE],
                                   wind_speed[i], wind_dir[i])
        damages[i] = compute_damage(result$burned, targets)
      }

      if (i %% 500 == 0) {
        cat(sprintf("  scenario %d / %d  (mean damage so far: %.2f)\n",
                    i, N, mean(damages[1:i])))
      }
    }

    results[[plan_name]] = data.frame(
      scenario   = seq_len(N),
      damage     = damages,
      wind_speed = wind_speed,
      wind_dir   = wind_dir,
      ign_row    = ignition[, "row"],
      ign_col    = ignition[, "col"]
    )

    cat(sprintf("  Done. Mean damage = %.2f\n", mean(damages)))
  }

  results
}


#' Compute risk measures from a damage vector
#'
#' Returns mean damage with 95% CI (normal approximation) and
#' CVaR at level alpha with 95% CI (bootstrap).
#'
#' @param damages Numeric vector of simulated damages.
#' @param alpha Numeric in (0,1). CVaR is mean of worst 1-alpha
#'   fraction. Default 0.9 (worst 10%).
#' @param n_boot Integer, bootstrap replicates for CVaR CI.
#' @return Data frame with columns measure, estimate, ci_lower, ci_upper.
compute_risk_measures = function(damages, alpha = 0.9, n_boot = 2000) {

  n = length(damages)

  ## Mean damage + 95% CI
  mu    = mean(damages)
  se_mu = sd(damages) / sqrt(n)
  ci_mu = mu + c(-1, 1) * qnorm(0.975) * se_mu

  ## CVaR_alpha = E[X | X >= VaR_alpha]
  cvar_point = function(x) {
    threshold = quantile(x, probs = alpha)
    mean(x[x >= threshold])
  }

  cvar_est = cvar_point(damages)

  ## bootstrap CI for CVaR
  boot_cvar = replicate(n_boot, {
    idx = sample.int(n, replace = TRUE)
    cvar_point(damages[idx])
  })
  ci_cvar = quantile(boot_cvar, probs = c(0.025, 0.975))

  data.frame(
    measure  = c("mean_damage", paste0("CVaR_", alpha)),
    estimate = c(mu, cvar_est),
    ci_lower = c(ci_mu[1], ci_cvar[[1]]),
    ci_upper = c(ci_mu[2], ci_cvar[[2]])
  )
}
