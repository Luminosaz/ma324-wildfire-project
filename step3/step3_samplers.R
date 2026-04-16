# ============================================================
# step3_samplers.R — Samplers for Monte Carlo fire evaluation
#
# Wind direction:  acceptance-rejection on von Mises mixture
# Wind speed:      Weibull(shape=2, scale=7) via rweibull
# Ignition:        uniform random vegetated cell
# ============================================================

#' Implements the wind-direction density f_Theta(theta) defined in the Model.
dwind_dir = function(theta) {
  p1  = 0.7;   mu1 = 1.5 * pi;  kappa1 = 3
  p2  = 0.3;   mu2 = 0.75 * pi; kappa2 = 2

  vm1 = exp(kappa1 * cos(theta - mu1)) / (2 * pi * besselI(kappa1, nu = 0))
  vm2 = exp(kappa2 * cos(theta - mu2)) / (2 * pi * besselI(kappa2, nu = 0))

  p1 * vm1 + p2 * vm2
}


#' Samples wind direction by acceptance-rejection with Uniform(0, 2*pi) proposal.
#' Envelope M = max f_Theta(theta) obtained by optimise().
#' Returns samples plus empirical and theoretical acceptance rates.
sample_wind_dir = function(n) {
  opt = optimise(dwind_dir, interval = c(0, 2 * pi), maximum = TRUE)
  M   = opt$objective

  samples      = numeric(n)
  n_accepted   = 0L
  n_total      = 0L
  total_accept = 0L

  while (n_accepted < n) {
    batch_size = max(n - n_accepted, 256L)

    theta = runif(batch_size, min = 0, max = 2 * pi)
    u     = runif(batch_size)

    accept       = u < dwind_dir(theta) / M
    n_total      = n_total + batch_size
    n_new        = sum(accept)
    total_accept = total_accept + n_new

    if (n_new == 0L) next

    n_take = min(n_new, n - n_accepted)
    samples[(n_accepted + 1):(n_accepted + n_take)] = theta[accept][1:n_take]
    n_accepted = n_accepted + n_take
  }

  list(
    samples          = samples,
    acceptance_rate  = total_accept / n_total,
    theoretical_rate = 1 / (2 * pi * M),
    M                = M
  )
}


#' Samples wind speed (km/h) from Weibull(shape = 2, scale = 7) by inverse transform.
sample_wind_speed = function(n) {
  rweibull(n, shape = 2, scale = 7)
}


#' Samples ignition cells uniformly over the vegetated cells of `landscape`.
sample_ignition = function(landscape, n) {
  veg_idx = which(landscape != 0, arr.ind = TRUE)
  chosen  = veg_idx[sample(nrow(veg_idx), size = n, replace = TRUE), , drop = FALSE]
  colnames(chosen) = c("row", "col")
  chosen
}
