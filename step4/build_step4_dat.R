# ============================================================
# build_step4_dat.R — Generate .dat for step4_reactive.mod
# ============================================================
# Calls edges_wind() for wind-dependent weights, then writes
# an AMPL .dat with: CELLS, IGNITION, TARGETS+weight,
# EDGES+w, PERM, R.
# ============================================================

source("edges_wind.R")


#' Build AMPL .dat file for the Step 4 reactive MILP
#'
#' Computes wind-dependent edge weights via \code{edges_wind()} and
#' writes a .dat matching \code{step4_reactive.mod}.
#'
#' @param landscape  Numeric matrix (21x21) from landscape.csv.
#' @param targets    Data frame with columns row, col, weight.
#' @param perm_cells Permanent cells from Step 3. Either an integer
#'   matrix (K x 2) with columns row, col, or a data frame with
#'   columns \code{row} and \code{col}.
#' @param wind_speed Scalar >= 0. Revealed wind speed (km/h).
#' @param wind_dir   Scalar in [0, 2*pi). Direction wind blows FROM.
#' @param alpha      Scalar in (0, 1]. Compression exponent for w.
#' @param R          Integer >= 0. Reactive clearing budget.
#' @param outfile    Character. Path for the .dat file.
#'
#' @return Invisible data frame of edges (with w column), for
#'   diagnostics.  Called for side effect of writing \code{outfile}.
build_step4_dat = function(landscape, targets, perm_cells,
                           wind_speed, wind_dir, alpha, R, outfile) {

  nr = nrow(landscape)
  nc = ncol(landscape)

  node_name = function(r, c) paste0(r, "_", c)

  # --- 1. set CELLS: all vegetated cells ----------------------
  cells = character(0)
  for (r in seq_len(nr)) {
    for (cc in seq_len(nc)) {
      if (landscape[r, cc] != 0) {
        cells = c(cells, node_name(r, cc))
      }
    }
  }

  # --- 2. set IGNITION: vegetated cells in row 21 -------------
  ignition = character(0)
  for (cc in seq_len(nc)) {
    if (landscape[nr, cc] != 0) {
      ignition = c(ignition, node_name(nr, cc))
    }
  }

  # --- 3. set TARGETS + param weight --------------------------
  target_nodes  = node_name(targets$row, targets$col)
  target_weight = targets$weight

  # --- 4. set PERM -------------------------------------------
  if (is.data.frame(perm_cells)) {
    perm_nodes = node_name(perm_cells$row, perm_cells$col)
  } else {
    perm_nodes = node_name(perm_cells[, 1], perm_cells[, 2])
  }

  # --- 5. set EDGES + param w (wind-dependent) ----------------
  edges = edges_wind(landscape, wind_speed = wind_speed,
                     wind_dir = wind_dir, alpha = alpha)

  edge_from = node_name(edges$from_row, edges$from_col)
  edge_to   = node_name(edges$to_row,   edges$to_col)
  keep      = (edge_from %in% cells) & (edge_to %in% cells)

  edge_from = edge_from[keep]
  edge_to   = edge_to[keep]
  w_vals    = edges$w[keep]

  # --- 6. Write .dat -----------------------------------------
  con = file(outfile, open = "w")
  on.exit(close(con))

  wl = function(...) writeLines(paste0(...), con)

  wl("# Auto-generated Step 4 reactive MILP data")
  wl("# wind_speed=", wind_speed,
     "  wind_dir=", round(wind_dir, 4),
     "  alpha=", alpha, "  R=", R)
  wl("")

  # set CELLS
  wl("set CELLS :=")
  wl("  ", paste(cells, collapse = " "))
  wl(";")
  wl("")

  # set IGNITION
  wl("set IGNITION :=")
  wl("  ", paste(ignition, collapse = " "))
  wl(";")
  wl("")

  # set TARGETS + param weight
  wl("param: TARGETS: weight :=")
  for (i in seq_along(target_nodes)) {
    wl("  ", target_nodes[i], "  ", target_weight[i])
  }
  wl(";")
  wl("")

  # set PERM
  wl("set PERM :=")
  if (length(perm_nodes) > 0) {
    wl("  ", paste(perm_nodes, collapse = " "))
  }
  wl(";")
  wl("")

  # set EDGES + param w
  wl("param: EDGES: w :=")
  for (i in seq_along(edge_from)) {
    wl("  ", edge_from[i], "  ", edge_to[i], "  ",
       format(w_vals[i], digits = 8))
  }
  wl(";")
  wl("")

  # param R
  wl("param R := ", R, ";")

  cat(sprintf("Wrote %s — %d cells, %d perm, %d edges, R=%d\n",
              outfile, length(cells), length(perm_nodes),
              length(edge_from), R))

  invisible(edges[keep, ])
}
