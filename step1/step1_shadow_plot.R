# =============================================================================
# Step 1: Shadow price aggregation and bottleneck heatmap
# =============================================================================

rm(list = ls())

landscape = as.matrix(read.table("../Project data and resources/landscape.csv", sep = ","))
targets   = read.csv("../Project data and resources/targets.csv")

# -----------------------------------------------------------------------------
# Aggregate shadow prices by receiving cell
# clearing cell j kills all edges into j, so score(j) = sum of shadow prices -> j
# -----------------------------------------------------------------------------

#' Aggregate shadow prices by receiving cell
#'
#' Reads a shadow price CSV and sums dual values on all edges pointing
#' into each cell. Clearing cell j removes all incoming edges, so the
#' sum represents how much max-flow would drop if j were cleared.
#'
#' @param csv_file path to shadow_prices CSV (from AMPL .run output)
#' @return 21x21 matrix of cell-level bottleneck scores
aggregate_shadow = function(csv_file) {
    df = read.csv(csv_file)

    # grid edges only, positive shadow price
    df = df[!df$from %in% c("S", "T") & !df$to %in% c("S", "T"), ]
    df = df[df$shadow_price > 1e-8, ]

    if (nrow(df) == 0) return(matrix(0, 21, 21))

    to_parts = strsplit(as.character(df$to), "_")
    df$to_row = as.integer(sapply(to_parts, `[`, 1))
    df$to_col = as.integer(sapply(to_parts, `[`, 2))

    agg = aggregate(shadow_price ~ to_row + to_col, data = df, FUN = sum)

    mat = matrix(0, 21, 21)
    for (k in seq_len(nrow(agg))) {
        mat[agg$to_row[k], agg$to_col[k]] = agg$shadow_price[k]
    }
    mat
}

# -----------------------------------------------------------------------------
# Heatmap panel
# -----------------------------------------------------------------------------

#' Plot a bottleneck heatmap for one target
#'
#' Draws a white-to-red heatmap of cell scores with blue outlines on targets.
#' Row 1 (north) is at the top.
#'
#' @param mat 21x21 matrix of bottleneck scores
#' @param title plot title (should state the finding)
#' @param target_cells data frame with row, col columns for blue outlines
plot_heatmap = function(mat, title, target_cells = NULL) {
    cols = colorRampPalette(c("white", "firebrick"))(256)
    max_val = max(mat)
    if (max_val < 1e-10) max_val = 1

    # flip so row 1 (north) is at top
    mat_flip = mat[21:1, ]

    image(x = 1:21, y = 1:21, z = t(mat_flip),
          col = cols, zlim = c(0, max_val), axes = FALSE,
          xlab = "Column", ylab = "Row", main = title)

    axis(1, at = seq(1, 21, by = 5))
    axis(2, at = seq(1, 21, by = 5), labels = seq(21, 1, by = -5))
    abline(h = 0:21 + 0.5, col = "grey85", lwd = 0.3)
    abline(v = 0:21 + 0.5, col = "grey85", lwd = 0.3)

    # target outlines in blue
    if (!is.null(target_cells)) {
        for (k in seq_len(nrow(target_cells))) {
            r = target_cells$row[k]
            cc = target_cells$col[k]
            y_flip = 22 - r
            rect(cc - 0.5, y_flip - 0.5, cc + 0.5, y_flip + 0.5,
                 border = "dodgerblue", lwd = 1.5)
        }
    }
    box()
}

# -----------------------------------------------------------------------------
# Load, aggregate, plot
# -----------------------------------------------------------------------------

mat_settle  = aggregate_shadow("shadow_prices_settlement.csv")
mat_wetland = aggregate_shadow("shadow_prices_wetland.csv")

settle_tgt  = targets[targets$weight == 10, ]
wetland_tgt = targets[targets$weight == 5, ]

png("figures/step1_shadow_heatmap.png", width = 2000, height = 1100, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
plot_heatmap(mat_settle,  "Settlement: pressure on cols 15-16", settle_tgt)
plot_heatmap(mat_wetland, "Wetland: pressure along cols 7-8 and row 7", wetland_tgt)
mtext("Fire pressure concentrates at target perimeters, not along current corridor",
      outer = TRUE, cex = 0.9, font = 2)
dev.off()

cat("Saved: figures/step1_shadow_heatmap.png\n")

# -----------------------------------------------------------------------------
# Top bottleneck cells
# -----------------------------------------------------------------------------

#' Print top n bottleneck cells to console
#'
#' @param mat 21x21 matrix of bottleneck scores
#' @param label name for display (e.g. "Settlement")
#' @param n how many cells to print (default 10)
print_top = function(mat, label, n = 10) {
    df = expand.grid(row = 1:21, col = 1:21)
    df$score = as.vector(mat)
    df = df[df$score > 1e-8, ]
    df = df[order(-df$score), ]
    cat(sprintf("\n--- Top %d bottleneck cells: %s ---\n", n, label))
    print(head(df, n), row.names = FALSE)
}

print_top(mat_settle,  "Settlement")
print_top(mat_wetland, "Wetland")
