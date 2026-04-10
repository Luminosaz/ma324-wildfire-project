# =============================================================================
# Step 1: Data exploration and edge list preparation
# =============================================================================

rm(list = ls())
source("../Project data and resources/fire-simulator.r")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

landscape = as.matrix(read.table("../Project data and resources/landscape.csv", sep = ","))
targets   = read.csv("../Project data and resources/targets.csv")

settlement = targets[targets$weight == 10, ]
wetland    = targets[targets$weight == 5, ]

cat("Landscape:", nrow(landscape), "x", ncol(landscape), "\n")
cat("Vegetated cells:", sum(landscape > 0), "\n")
cat("Bare ground:", sum(landscape == 0), "\n")

# Quick sanity check: settlement NE, wetland NW
cat("Settlement rows:", range(settlement$row), " cols:", range(settlement$col), "\n")
cat("Wetland rows:",    range(wetland$row),    " cols:", range(wetland$col), "\n")

# -----------------------------------------------------------------------------
# Landscape overview plot
# -----------------------------------------------------------------------------

# Current corridor: row 10, cols 3-19 (17 cells, 2 open on each side)
corridor = cbind(rep(10, 17), 3:19)

png("figures/landscape_overview.png", width = 1200, height = 1200, res = 200)
plot_fire(landscape, firebreaks = corridor, targets = targets,
          main = "Landscape overview: targets and current corridor")
text(10, 21 - 4.5 + 0.5, "Wetland\n(w=5)", cex = 0.6, font = 2, col = "white")
text(17.5, 21 - 4.5 + 0.5, "Settlement\n(w=10)", cex = 0.6, font = 2, col = "white")
text(11, 21 - 10 + 0.5, "Current corridor (row 10)", cex = 0.5, col = "black", pos = 3)
dev.off()

# -----------------------------------------------------------------------------
# Spread probabilities (no wind, for Step 1 max-flow)
# -----------------------------------------------------------------------------

edges = spread_probabilities(landscape, wind_speed = 0, wind_dir = 0)
cat("\nEdges:", nrow(edges), "\n")
cat("p_burn range:", round(range(edges$p_burn), 3), "\n")
print(round(quantile(edges$p_burn, c(0, 0.25, 0.5, 0.75, 1)), 3))

saveRDS(edges, "edges_nowind.rds")
