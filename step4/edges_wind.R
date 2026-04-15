# ============================================================
# edges_wind.R — Wind-dependent edge weights for reactive MILP
# ============================================================
# Wraps the simulator's spread_probabilities() to produce an
# edge data frame augmented with compressed weight w = p_burn^alpha.
# Used by Step 4 to generate wind-specific .dat files for the
# reactive firebreak MILP.
# ============================================================


#' Compute wind-dependent edge weights for the MILP
#'
#' Calls \code{spread_probabilities()} from fire-simulator.r and
#' appends a compressed weight column \code{w = p_burn^alpha} suitable
#' for the continuous-reachability MILP (Model B).
#'
#' @param landscape  Numeric matrix (21x21). Cell encoding: tens digit =
#'   vegetation type (1--3), units digit = density (1--3), 0 = bare ground.
#' @param wind_speed Scalar >= 0. Wind speed in km/h.
#' @param wind_dir   Scalar in [0, 2*pi). Direction the wind blows FROM,
#'   following the simulator convention. Do NOT add pi — the caller is
#'   responsible for passing the raw meteorological direction.
#' @param alpha      Scalar in (0, 1]. Compression exponent. Smaller values
#'   compress the range of edge weights, preventing products from decaying
#'   below solver tolerance over many hops. No default — caller must specify.
#'
#' @return A data frame identical to the output of
#'   \code{spread_probabilities()} (columns \code{from_row}, \code{from_col},
#'   \code{to_row}, \code{to_col}, \code{p_burn}) plus an additional column
#'   \code{w = p_burn^alpha}.
#'
#' @details
#'   The function assumes \code{spread_probabilities()} is already loaded
#'   (sourced from \code{fire-simulator.r}). It performs no wind-direction
#'   transformation — the value is passed straight through.
#'
#' @examples
#' \dontrun{
#'   source("fire-simulator.r")
#'   landscape = as.matrix(read.csv("landscape.csv", header = FALSE))
#'   edges = edges_wind(landscape, wind_speed = 10, wind_dir = 3*pi/2,
#'                      alpha = 0.5)
#'   head(edges)
#' }
edges_wind = function(landscape, wind_speed, wind_dir, alpha) {

  # --- input checks --------------------------------------------------
  stopifnot(is.matrix(landscape), is.numeric(landscape))
  stopifnot(length(wind_speed) == 1L, wind_speed >= 0)
  stopifnot(length(wind_dir)   == 1L, wind_dir >= 0, wind_dir < 2 * pi)
  stopifnot(length(alpha)      == 1L, alpha > 0, alpha <= 1)

  # --- call the simulator's edge builder -----------------------------
  edges = spread_probabilities(landscape,
                               wind_speed = wind_speed,
                               wind_dir   = wind_dir)

  # --- append compressed weight --------------------------------------
  edges$w = edges$p_burn ^ alpha

  edges
}
