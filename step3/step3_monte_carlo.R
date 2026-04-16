# ============================================================
# step3_monte_carlo.R — Monte Carlo evaluation of firebreak
#   plans under uncertain wind and ignition.
#
# Dependencies:
#   source("step3_samplers.R")
#   source("../../Project data and resources/fire-simulator.r")
# ============================================================


#' Evaluates all plans on a pre-sampled batch of (wind_dir, wind_speed, ignition)
#' scenarios. Accepting the scenarios as arguments lets the caller reuse the
#' same draws for diagnostics and the Monte Carlo evaluation.
#' Ignition on a cleared cell yields zero damage (no fire starts).
run_mc_evaluation = function(landscape, targets, plans,
                             wind_dir, wind_speed, ignition) {

  N = length(wind_dir)
  stopifnot(length(wind_speed) == N, nrow(ignition) == N)

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


#' Computes mean damage, a 95% t-interval, and bootstrap CI for CVaR_alpha.
compute_risk_measures = function(damages, alpha = 0.9, n_boot = 2000) {

  n = length(damages)

  ## Mean damage + 95% t-interval
  mu    = mean(damages)
  se_mu = sd(damages) / sqrt(n)
  ci_mu = mu + c(-1, 1) * qt(0.975, df = n - 1) * se_mu

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
