# =============================================================================
# Step 2: Generate AMPL .dat files for firebreak MILP models
# =============================================================================

rm(list = ls())
source("../../Project data and resources/fire-simulator.r")

#' Generate AMPL .dat file for Step 2 firebreak MILP
#'
#' @param edges_df  Data frame from spread_probabilities() with columns:
#'                  from_row, from_col, to_row, to_col, p_burn
#' @param targets   Data frame from targets.csv with columns: row, col, weight
#' @param outfile   Path to write the .dat file
#' @param landscape Integer matrix (21x21) of vegetation codes
#'
#' @return Invisible NULL. Called for side effect of writing the .dat file.
generate_step2_dat <- function(edges_df, targets, outfile, landscape) {

  node_name <- function(r, c) paste0(r, "_", c)

  nr <- nrow(landscape)
  nc <- ncol(landscape)

  # --- 1. set CELLS: all vegetated (non-bare) cells ----
  cells <- c()
  for (r in 1:nr) {
    for (c in 1:nc) {
      if (landscape[r, c] != 0) {
        cells <- c(cells, node_name(r, c))
      }
    }
  }

  # --- 2. set IGNITION: vegetated cells in row 21 ----
  ignition <- c()
  for (c in 1:nc) {
    if (landscape[nr, c] != 0) {
      ignition <- c(ignition, node_name(nr, c))
    }
  }

  # --- 3. set TARGETS + param weight ----
  target_nodes  <- node_name(targets$row, targets$col)
  target_weight <- targets$weight

  # --- 4. set EDGES + param p_burn ----
  edge_from <- node_name(edges_df$from_row, edges_df$from_col)
  edge_to   <- node_name(edges_df$to_row,   edges_df$to_col)
  keep      <- (edge_from %in% cells) & (edge_to %in% cells)

  edge_from <- edge_from[keep]
  edge_to   <- edge_to[keep]
  p_burn    <- edges_df$p_burn[keep]

  # --- 5. Write .dat file ----
  con <- file(outfile, open = "w")
  on.exit(close(con))

  wl <- function(...) writeLines(paste0(...), con)

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

  # set EDGES + param p_burn
  wl("param: EDGES: p_burn :=")
  for (i in seq_along(edge_from)) {
    wl("  ", edge_from[i], "  ", edge_to[i], "  ", format(p_burn[i], digits = 8))
  }
  wl(";")
  wl("")

  cat(sprintf("Wrote %s: %d cells, %d ignition, %d targets, %d edges\n",
              outfile, length(cells), length(ignition),
              length(target_nodes), length(edge_from)))
}

# -----------------------------------------------------------------------------
# Load data and generate .dat
# -----------------------------------------------------------------------------

landscape <- as.matrix(read.table(
  "../../Project data and resources/landscape.csv", sep = ","
))
targets <- read.csv("../../Project data and resources/targets.csv")
edges   <- spread_probabilities(landscape, wind_speed = 0, wind_dir = 0)

generate_step2_dat(edges, targets, outfile = "step2.dat",
                   landscape = landscape)
