# =============================================================================
# Step 1: Generate AMPL .dat files for max-flow model
# =============================================================================

rm(list = ls())
source("../Project data and resources/fire-simulator.r")

# -----------------------------------------------------------------------------
# Generate .dat file for maxflow.mod
# -----------------------------------------------------------------------------

generate_maxflow_dat = function(edges_df, target_cells, outfile,
                                ignition_row = 21, big_M = 999) {
    node_name = function(r, c) paste0(r, "_", c)

    # All grid nodes from edge list (bare ground already excluded by spread_probabilities)
    grid_nodes = unique(c(
        node_name(edges_df$from_row, edges_df$from_col),
        node_name(edges_df$to_row,   edges_df$to_col)
    ))
    all_nodes = c("S", "T", sort(grid_nodes))

    # Ignition: all cells in ignition row that exist as nodes
    ignition_nodes = grep(paste0("^", ignition_row, "_"), grid_nodes, value = TRUE)
    target_nodes   = intersect(node_name(target_cells$row, target_cells$col), grid_nodes)

    if (length(ignition_nodes) == 0) stop("No ignition nodes found")
    if (length(target_nodes) == 0)   stop("No target nodes found")

    # Build arc list
    arcs = data.frame(
        from = node_name(edges_df$from_row, edges_df$from_col),
        to   = node_name(edges_df$to_row,   edges_df$to_col),
        cap  = edges_df$p_burn
    )

    source_arcs = data.frame(from = "S", to = ignition_nodes, cap = big_M)
    sink_arcs   = data.frame(from = target_nodes, to = "T",   cap = big_M)
    return_arc  = data.frame(from = "T", to = "S", cap = big_M)

    all_arcs = rbind(return_arc, source_arcs, sink_arcs, arcs)

    # Write .dat
    f = file(outfile, open = "w")
    on.exit(close(f))

    writeLines(sprintf("# Max-flow data | ignition row: %d | targets: %d | arcs: %d",
                       ignition_row, length(target_nodes), nrow(arcs)), f)
    writeLines("", f)

    # Nodes
    writeLines("set NODES :=", f)
    chunks = split(all_nodes, ceiling(seq_along(all_nodes) / 10))
    for (chunk in chunks) {
        writeLines(paste("   ", paste(chunk, collapse = "  ")), f)
    }
    writeLines(";", f)
    writeLines("", f)

    # Arcs + capacities
    writeLines("param: ARCS: cap :=", f)
    for (k in seq_len(nrow(all_arcs))) {
        writeLines(sprintf("    %-8s %-8s  %.6f", all_arcs$from[k], all_arcs$to[k], all_arcs$cap[k]), f)
    }
    writeLines(";", f)

    cat(sprintf("Written %s: %d nodes, %d arcs\n", outfile, length(all_nodes), nrow(all_arcs)))
}

# -----------------------------------------------------------------------------
# Load data and generate both .dat files
# -----------------------------------------------------------------------------

landscape = as.matrix(read.table("../Project data and resources/landscape.csv", sep = ","))
targets   = read.csv("../Project data and resources/targets.csv")
edges     = spread_probabilities(landscape, wind_speed = 0, wind_dir = 0)

settlement = targets[targets$weight == 10, ]
wetland    = targets[targets$weight == 5, ]

generate_maxflow_dat(edges, settlement, "step1_settlement.dat")
generate_maxflow_dat(edges, wetland,    "step1_wetland.dat")
