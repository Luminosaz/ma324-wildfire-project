# ============================================================
# test_edges_wind.R — Sanity checks for edges_wind()
# ============================================================
# Two checks:
#   1. V = 0  : structural integrity (p_burn, w in [0,1],
#               no self-loops, no duplicated directed pairs)
#   2. V = 10, wind_dir = pi/2 : from cell (10,10) the
#               westward neighbour gets higher p_burn than
#               the eastward one
#
# Requires: fire-simulator.r, landscape.csv (in ../../Project data
#           and resources/), edges_wind.R (in working directory)
# Runs from: code/step4/
# ============================================================

cat("=== test_edges_wind.R ===\n\n")

# --- setup -------------------------------------------------------
# Runs from code/step4/; simulator and data live two levels up
data_dir = "../../Project data and resources"
source(file.path(data_dir, "fire-simulator.r"))
source("edges_wind.R")

landscape = as.matrix(read.csv(file.path(data_dir, "landscape.csv"),
                               header = FALSE))
n_pass = 0L
n_fail = 0L

check = function(desc, expr) {
  ok = tryCatch(expr, error = function(e) FALSE)
  if (isTRUE(ok)) {
    cat(sprintf("  PASS: %s\n", desc))
    n_pass <<- n_pass + 1L
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    n_fail <<- n_fail + 1L
  }
}

# =================================================================
# CHECK 1 — structural integrity at V = 0 (no wind)
# =================================================================
cat("Check 1: V = 0, no wind\n")

e0 = edges_wind(landscape, wind_speed = 0, wind_dir = 0, alpha = 0.5)

# 1a. p_burn in [0, 1]
check("all p_burn in [0, 1]",
      all(e0$p_burn >= 0 & e0$p_burn <= 1))

# 1b. w in [0, 1]
check("all w in [0, 1]",
      all(e0$w >= 0 & e0$w <= 1))

# 1c. w == p_burn^alpha
check("w equals p_burn^0.5 (to machine precision)",
      all(abs(e0$w - e0$p_burn^0.5) < 1e-12))

# 1d. no self-loops (from_row == to_row AND from_col == to_col)
self_loop = (e0$from_row == e0$to_row) & (e0$from_col == e0$to_col)
check("no self-loops",
      !any(self_loop))

# 1e. no duplicated directed pairs
edge_key = paste(e0$from_row, e0$from_col, e0$to_row, e0$to_col, sep = "_")
check("no duplicated directed pairs",
      !any(duplicated(edge_key)))

cat("\n")

# =================================================================
# CHECK 2 — wind directionality: wind_dir = pi/2, V = 10
# =================================================================
cat("Check 2: V = 10, wind_dir = pi/2\n")

ew = edges_wind(landscape, wind_speed = 10, wind_dir = pi / 2, alpha = 0.5)

# Find edges FROM cell (10, 10) to its east and west neighbours
from_10_10 = (ew$from_row == 10) & (ew$from_col == 10)

# Westward neighbour: (10, 9)
west = from_10_10 & (ew$to_row == 10) & (ew$to_col == 9)
# Eastward neighbour: (10, 11)
east = from_10_10 & (ew$to_row == 10) & (ew$to_col == 11)

if (sum(west) != 1L || sum(east) != 1L) {
  cat("  FAIL: could not find unique west/east edges from (10,10)\n")
  cat(sprintf("         west edges found: %d, east edges found: %d\n",
              sum(west), sum(east)))
  n_fail = n_fail + 1L
} else {
  p_west = ew$p_burn[west]
  p_east = ew$p_burn[east]
  cat(sprintf("  p_burn to west (10,9)  = %.6f\n", p_west))
  cat(sprintf("  p_burn to east (10,11) = %.6f\n", p_east))

  check("westward p_burn > eastward p_burn",
        p_west > p_east)
}

cat("\n")

# =================================================================
# Summary
# =================================================================
cat(sprintf("Results: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0L) {
  stop("Some checks failed — see above.")
} else {
  cat("All checks passed.\n")
}
