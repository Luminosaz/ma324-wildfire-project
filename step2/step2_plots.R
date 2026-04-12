# ============================================================
# step2_plots.R — Multi-panel comparison figure for Step 2
# ============================================================
# Produces a 3-panel PNG (200 DPI) for the report:
#   Panel 1: Model B damage vs budget at different α
#   Panel 2: Model A (τ=0) vs Model B (α=0.1) vs simulator (log scale)
#   Panel 3: Firebreak location maps for corridor, A B=16, B B=16
# ============================================================

source("../.Rprofile")

# --- Load data -----------------------------------------------
res1 = read.csv("step2_results.csv",   stringsAsFactors = FALSE)
res2 = read.csv("step2_results_A.csv", stringsAsFactors = FALSE)
sim = read.csv("step2_sim_results_K1000.csv", stringsAsFactors = FALSE)

data_dir  ="../../Project data and resources"
landscape = as.matrix(read.table(file.path(data_dir, "landscape.csv"),
                                  sep = ",", header = FALSE))
targets  = read.csv(file.path(data_dir, "targets.csv"))

all_res = rbind(res1, res2)

# =============================================================
# Helpers
# =============================================================

#' Parse semicolon-separated "r_c" strings into (row, col) matrix
#'
#' @param cell_str Single string of semicolon-separated node names
#'
#' @return Integer matrix with columns named `row` and `col`,
#'   or a 0-row matrix if input is NA or empty
parse_cells = function(cell_str) {
  if (is.na(cell_str) || nchar(trimws(cell_str)) == 0) {
    return(matrix(integer(0), ncol = 2,
                  dimnames = list(NULL, c("row", "col"))))
  }
  tokens = strsplit(cell_str, ";")[[1]]
  parts = strsplit(tokens, "_")
  rc    = do.call(rbind, lapply(parts, as.integer))
  colnames(rc) = c("row", "col")
  rc
}

#' Draw a 21x21 grid heatmap showing landscape, targets, and firebreaks
#'
#' @param landscape  21x21 integer matrix (vegetation encoding)
#' @param targets    Data frame with row, col, weight columns
#' @param breaks_rc  Two-column (row, col) matrix of cleared cells, or NULL
#' @param title      Character string for plot title
#'
#' @return NULL (called for side effect of drawing)
draw_grid = function(landscape, targets, breaks_rc, title) {
  nr = nrow(landscape)
  nc = ncol(landscape)

  # Vegetation type: tens digit (1=grass, 2=shrub, 3=forest, 0=bare)
  veg = landscape %/% 10
  # Colour map: bare=white, grass=lightyellow, shrub=khaki, forest=darkgreen
  veg_cols = c("0" = "grey95", "1" = "#d4e89c", "2" = "#8db944", "3" = "#2d6a1e")

  # Build colour matrix
  col_mat = matrix("grey95", nr, nc)
  for (r in 1:nr) {
    for (c in 1:nc) {
      v = as.character(veg[r, c])
      if (v %in% names(veg_cols)) col_mat[r, c] = veg_cols[v]
    }
  }

  # Mark targets
  for (i in seq_len(nrow(targets))) {
    col_mat[targets$row[i], targets$col[i]] <-
      if (targets$weight[i] == 10) "#4a90d9" else "#9b59b6"
  }

  # Mark firebreaks
  if (!is.null(breaks_rc) && nrow(breaks_rc) > 0) {
    for (i in seq_len(nrow(breaks_rc))) {
      col_mat[breaks_rc[i, 1], breaks_rc[i, 2]] ="#e74c3c"
    }
  }

  # Plot as image (row 1 = top → need to flip for image())
  plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5),
       xlab = "col", ylab = "row", asp = 1, main = title,
       cex.main = 0.9, xaxt = "n", yaxt = "n")
  axis(1, at = seq(1, nc, 4), cex.axis = 0.7)
  axis(2, at = seq(1, nr, 4), cex.axis = 0.7, las = 1)

  for (r in 1:nr) {
    for (c in 1:nc) {
      rect(c - 0.5, r - 0.5, c + 0.5, r + 0.5,
           col = col_mat[r, c], border = NA)
    }
  }
  # Light grid
  abline(h = 1:(nr) - 0.5, col = "grey80", lwd = 0.3)
  abline(v = 1:(nc) - 0.5, col = "grey80", lwd = 0.3)
  box()
}

# =============================================================
# Panel 1 data: Model B damage vs B at each α
# =============================================================
modelB = all_res[all_res$model == "B", ]
alpha_vals = sort(unique(modelB$param_value))
B_vals     = sort(unique(modelB$B))

# Colours for α lines
alpha_cols = c("#e74c3c", "#e67e22", "#2ecc71", "#3498db", "#9b59b6")
names(alpha_cols) = as.character(alpha_vals)

# =============================================================
# Panel 2 data: Model A τ=0 vs Model B α=0.1 vs simulator
# =============================================================
modelA_tau0  = all_res[all_res$model == "A" &
                        all_res$param_value == 0, ]
modelA_tau0  = modelA_tau0[!duplicated(modelA_tau0$B), ]
modelA_tau0  = modelA_tau0[order(modelA_tau0$B), ]

modelB_a01   = modelB[modelB$param_value == 0.1, ]
modelB_a01   = modelB_a01[!duplicated(modelB_a01$B), ]
modelB_a01   = modelB_a01[order(modelB_a01$B), ]

# Simulator results (exclude corridor and baseline for overlay)
sim_milp = sim[sim$model %in% c("A", "B"), ]

# =============================================================
# Panel 3 data: key plans — use B=16 where A and B disagree
# =============================================================
# Model A B=16 τ=0 (uses only 13 of 16 cells — settlement only)
a16 = all_res[all_res$model == "A" & all_res$B == 16 &
               all_res$param_value == 0 &
               nchar(trimws(all_res$cleared_cells)) > 0, ]
a16_rc = if (nrow(a16) > 0) parse_cells(a16$cleared_cells[1]) else NULL

# Model B B=16 α=0.1 (uses all 16 cells — adds wetland protection)
b16 = modelB[modelB$B == 16 & modelB$param_value == 0.1 &
              nchar(trimws(modelB$cleared_cells)) > 0, ]
b16_rc = if (nrow(b16) > 0) parse_cells(b16$cleared_cells[1]) else NULL

# Current corridor: row 10, cols 3–19
corr_rc = cbind(row = rep(10, 17), col = 3:19)

# =============================================================
# Draw figure
# =============================================================
png("step2_comparison.png", width = 11, height = 8.5,
    units = "in", res = 200)

layout(matrix(c(1, 1, 1, 2, 2, 2,
                3, 3, 4, 4, 5, 5), nrow = 2, byrow = TRUE),
       heights = c(1, 1.2))
par(cex.lab = 0.95, cex.axis = 0.85)

# --- Panel 1: Model B damage vs B at different α -------------
par(mar = c(4, 4.2, 2.5, 1))
plot(NA, xlim = range(B_vals), ylim = c(0, max(modelB$damage) * 1.05),
     xlab = "Budget B", ylab = "MILP damage",
     main = "Model B: damage vs budget")
for (i in seq_along(alpha_vals)) {
  a  = alpha_vals[i]
  dd = modelB[modelB$param_value == a, ]
  dd = dd[order(dd$B), ]
  lines(dd$B, dd$damage, col = alpha_cols[i], lwd = 2, type = "b",
        pch = 15 + i - 1, cex = 0.9)
}
legend("topright", legend = paste0("\u03b1=", alpha_vals),
       col = alpha_cols, lwd = 2, pch = 15:19,
       cex = 0.75, bg = "white", ncol = 1)

# --- Panel 2: MILP predictions + simulator on same scale ------
# Use log-ish y-axis: plot simulator range properly,
# show MILP as inset or use two sub-axes.
# Approach: simulator-scale plot with MILP on secondary axis.
par(mar = c(4, 4.2, 2.5, 4.2))

# Simulator scale
sim_ymax = max(sim$sim_ci_hi[!is.na(sim$sim_ci_hi)], na.rm = TRUE) * 1.3
plot(NA, xlim = range(B_vals), ylim = c(0, sim_ymax),
     xlab = "Budget B", ylab = "Simulated damage (K=1000)",
     main = "MILP prediction vs simulator")

# Baseline and corridor reference lines
bl_mean  = sim$sim_mean[sim$plan_label == "baseline"]
corr_mean = sim$sim_mean[sim$plan_label == "corridor_row10"]
abline(h = bl_mean,   lty = 3, col = "grey50", lwd = 1.2)
abline(h = corr_mean, lty = 2, col = "grey50", lwd = 1.2)
text(19, bl_mean + 0.7,   "no firebreaks", cex = 0.65, col = "grey40", adj = 1)
text(19, corr_mean - 0.7, "corridor",       cex = 0.65, col = "grey40", adj = 1)

# Simulator points with error bars
for (i in seq_len(nrow(sim_milp))) {
  s   = sim_milp[i, ]
  pcol = if (s$model == "A") "#c0392b" else "#e67e22"
  ppch = if (s$model == "A") 17 else 15
  # Slight x-jitter to separate A and B at same budget
  xj   = s$B + if (s$model == "A") -0.25 else 0.25
  points(xj, s$sim_mean, pch = ppch, col = pcol, cex = 1.1, lwd = 1.5)
  arrows(xj, s$sim_ci_lo, xj, s$sim_ci_hi,
         code = 3, angle = 90, length = 0.04, col = pcol, lwd = 1.3)
}

# MILP predictions on secondary y-axis (right)
milp_max = max(c(modelA_tau0$damage, modelB_a01$damage), na.rm = TRUE)
scale_f = sim_ymax / milp_max
lines(modelA_tau0$B, modelA_tau0$damage * scale_f,
      col = "#2c3e50", lwd = 1.5, lty = 2)
points(modelA_tau0$B, modelA_tau0$damage * scale_f,
       pch = 2, col = "#2c3e50", cex = 0.8)
lines(modelB_a01$B, modelB_a01$damage * scale_f,
      col = "#2980b9", lwd = 1.5, lty = 2)
points(modelB_a01$B, modelB_a01$damage * scale_f,
       pch = 0, col = "#2980b9", cex = 0.8)

axis(4, at = pretty(c(0, sim_ymax), 5),
     labels = round(pretty(c(0, sim_ymax), 5) / scale_f),
     cex.axis = 0.75, col.axis = "grey30")
mtext("MILP damage (dashed)", side = 4, line = 2.8,
      cex = 0.7, col = "grey30")

legend("topright",
       legend = c("A sim \u00b1 95% CI", "B sim \u00b1 95% CI",
                  expression("A " * tau * "=0 (MILP)"),
                  expression("B " * alpha * "=0.1 (MILP)")),
       col = c("#c0392b", "#e67e22", "#2c3e50", "#2980b9"),
       pch = c(17, 15, 2, 0), lty = c(NA, NA, 2, 2),
       lwd = c(NA, NA, 1.5, 1.5),
       cex = 0.65, bg = "white")

# --- Panels 3-5: firebreak location maps ---------------------
par(mar = c(3, 3.5, 2.5, 1))

draw_grid(landscape, targets, corr_rc,
          "Corridor (row 10, 17 cells)")
draw_grid(landscape, targets, a16_rc,
          expression("Model A  B=16  " * tau * "=0  (13 used)"))
draw_grid(landscape, targets, b16_rc,
          expression("Model B  B=16  " * alpha * "=0.1  (16 used)"))

# Shared legend for grid maps
par(xpd = NA)
legend("bottom", inset = c(0, -0.18), horiz = TRUE,
       legend = c("Grass", "Shrub", "Forest",
                  "Settlement", "Wetland", "Firebreak"),
       fill = c("#d4e89c", "#8db944", "#2d6a1e",
                "#4a90d9", "#9b59b6", "#e74c3c"),
       cex = 0.6, border = NA, bty = "n")

dev.off()
cat("Saved step2_comparison.png\n")
