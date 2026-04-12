# MA324 Report Draft

Working draft for Cadmus. This will keep updating as we progress through the steps.

---

## Current Solution

(to write last — needs full picture from all steps)

---

## Proposal

(to write last — after all steps done)

---

## Scaffolding

### Step 1: Network Bottleneck Model

**Model**

For Step 1, we model the fire grid as a network flow problem. Each vegetated cell on the 21×21 grid becomes a node, and we draw directed edges between adjacent cells (including diagonals, so 8 directions). The capacity of each edge is set to the spread probability p_burn from the simulator with no wind.

We add a source node S that connects to all cells in row 21 (the south edge where fire starts), and a sink node T that connects from the target cells. There is also a return arc from T back to S so we can use the standard max-flow formulation.

The LP maximises total flow from S to T:

$$\begin{aligned} \max \quad & \text{flow}(T, S) \\ \text{s.t.} \quad & \sum_{(i,v) \in E} \text{flow}(i,v) = \sum_{(v,j) \in E} \text{flow}(v,j) \quad \forall \, v \in V \\ & \text{flow}(i,j) \leq \text{cap}(i,j) \quad \forall \, (i,j) \in E \end{aligned}$$

We write the capacity constraint as a named constraint in AMPL so we can read shadow prices with .dual. We run this twice — once for settlement as the target, once for wetland.

To figure out which cells are most important to clear, we look at the shadow prices. Since clearing a cell removes all edges going into it, we sum up the shadow prices of all incoming edges for each cell. A higher sum means that cell is a bigger bottleneck.

**Code**

Max-flow on fire spread network:

```ampl
# Sets
set NODES;
set ARCS within {NODES, NODES};

# Parameters
param cap {ARCS} >= 0;

# Variables — flow on each arc, lower bound only
var flow {(i,j) in ARCS} >= 0;

# Objective: maximise flow on return arc
maximize TotalFlow: flow["T","S"];

# Flow conservation at every node
subject to FlowBalance {v in NODES}:
    sum {(i,v) in ARCS} flow[i,v] - sum {(v,j) in ARCS} flow[v,j] = 0;

# Capacity as named constraint (not variable bound) for shadow prices
subject to Capacity {(i,j) in ARCS: i <> "T" or j <> "S"}:
    flow[i,j] <= cap[i,j];
```

R function to generate the AMPL .dat file from the spread probability edge list:

```r
#' Generate AMPL .dat file for the max-flow model
#'
#' Builds the network (source -> ignition row, spread edges, targets -> sink)
#' and writes it in AMPL-readable format matching step1_maxflow.mod.
#'
#' @param edges_df data frame from spread_probabilities() with columns
#'   from_row, from_col, to_row, to_col, p_burn
#' @param target_cells data frame with row, col columns (settlement or wetland)
#' @param outfile path for the output .dat file
#' @param ignition_row which row is the fire source (default 21 = south)
#' @param big_M large capacity for source/sink arcs (default 999)
#' @return writes .dat file, prints summary to console
generate_maxflow_dat = function(edges_df, target_cells, outfile,
                                ignition_row = 21, big_M = 999) {
    node_name = function(r, c) paste0(r, "_", c)

    grid_nodes = unique(c(
        node_name(edges_df$from_row, edges_df$from_col),
        node_name(edges_df$to_row,   edges_df$to_col)
    ))
    all_nodes = c("S", "T", sort(grid_nodes))

    ignition_nodes = grep(paste0("^", ignition_row, "_"), grid_nodes, value = TRUE)
    target_nodes   = intersect(node_name(target_cells$row, target_cells$col), grid_nodes)

    if (length(ignition_nodes) == 0) stop("No ignition nodes found")
    if (length(target_nodes) == 0)   stop("No target nodes found")

    arcs = data.frame(
        from = node_name(edges_df$from_row, edges_df$from_col),
        to   = node_name(edges_df$to_row,   edges_df$to_col),
        cap  = edges_df$p_burn
    )

    source_arcs = data.frame(from = "S", to = ignition_nodes, cap = big_M)
    sink_arcs   = data.frame(from = target_nodes, to = "T",   cap = big_M)
    return_arc  = data.frame(from = "T", to = "S", cap = big_M)

    all_arcs = rbind(return_arc, source_arcs, sink_arcs, arcs)

    f = file(outfile, open = "w")
    on.exit(close(f))

    writeLines(sprintf("# Max-flow data | ignition row: %d | targets: %d | arcs: %d",
                       ignition_row, length(target_nodes), nrow(arcs)), f)
    writeLines("", f)

    writeLines("set NODES :=", f)
    chunks = split(all_nodes, ceiling(seq_along(all_nodes) / 10))
    for (chunk in chunks) {
        writeLines(paste("   ", paste(chunk, collapse = "  ")), f)
    }
    writeLines(";", f)
    writeLines("", f)

    writeLines("param: ARCS: cap :=", f)
    for (k in seq_len(nrow(all_arcs))) {
        writeLines(sprintf("    %-8s %-8s  %.6f", all_arcs$from[k], all_arcs$to[k], all_arcs$cap[k]), f)
    }
    writeLines(";", f)

    cat(sprintf("Written %s: %d nodes, %d arcs\n", outfile, length(all_nodes), nrow(all_arcs)))
}
```

**Results**

We solved the max-flow for both targets:

| Target | Max flow |
|--------|----------|
| Settlement | 4.06 |
| Wetland | 3.795 |

Settlement has higher max-flow, meaning more fire can reach it — it is harder to protect.

[INSERT: figures/step1_shadow_heatmap.png]

**Figure 1:** Bottleneck pressure by cell for settlement (left) and wetland (right). Red intensity shows the sum of shadow prices on incoming edges. Settlement bottlenecks cluster around columns 15-16 and row 7, while wetland pressure concentrates at columns 7-8. Both are near the target perimeters, not along the current row 10 corridor.

**Discussion**

The current east-west corridor along row 10 does not coincide with the bottleneck corridors. The shadow prices show fire pressure is highest at rows 6-7, right next to the targets. The corridor is too far south and leaves two cells open on each side, so fire can bypass it.

Settlement is harder to protect than wetland, with a higher max-flow of 4.06 compared to 3.795. This is because the east side of the grid has denser forest vegetation, which gives higher spread probabilities toward the settlement.

This model treats spread probabilities as fixed capacities, so it ignores the stochastic nature of fire spread. It also assumes no wind and a fixed southern ignition front. In practice, wind can shift the direction of fire flow significantly, and ignition can occur anywhere. These limitations are addressed in Steps 2 and 3.

Based on the shadow prices, we would propose clearing cells along column 15 (rows 3-7) to protect the settlement, and along column 7 (rows 3-7) for the wetland. These form barriers at the points where fire pressure is highest. However, this does not account for budget constraints, which we introduce in Step 2.

---

### Step 2: Budget-Constrained Firebreak Optimisation

(todo)

---

### Step 3: Evaluation Under Uncertain Wind

(todo)

---

### Step 4: Reactive Response

(todo)

---

## Discussion

(to write after all steps — synthesise across steps)

---

## AI Declaration

(to write last)
