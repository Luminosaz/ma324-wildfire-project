#' Run K fire simulations and return mean damage with 95% CI
#'
#' @param landscape  21x21 matrix (firebreaks already applied as 0)
#' @param ignition   Two-column matrix (row, col) of ignition cells
#' @param targets    Data frame with row, col, weight columns
#' @param K          Number of replications
#'
#' @return Named list: mean_damage, se, ci_lo, ci_hi
run_simulation = function(landscape, ignition, targets, K) {
  damages = numeric(K)
  for (k in seq_len(K)) {
    result    = simulate_fire(landscape, ignition,
                                wind_speed = 0, wind_dir = 0)
    damages[k] = compute_damage(result$burned, targets)
  }
  mu = mean(damages)
  se = sd(damages) / sqrt(K)
  list(mean_damage = mu,
       se         = se,
       ci_lo      = mu - 1.96 * se,
       ci_hi      = mu + 1.96 * se)
}
