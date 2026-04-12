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

The .run script solves both targets with Gurobi and exports shadow prices to CSV via `Capacity[i,j].dual`.

R function to build the network and generate .dat files:

```r
#' Generate AMPL .dat file for the max-flow model
#'
#' @param edges_df data frame from spread_probabilities()
#' @param target_cells data frame with row, col columns
#' @param outfile path for the output .dat file
#' @param ignition_row fire source row (default 21)
#' @param big_M capacity for source/sink arcs (default 999)
#' @return writes .dat file to disk
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

    arcs = data.frame(
        from = node_name(edges_df$from_row, edges_df$from_col),
        to   = node_name(edges_df$to_row,   edges_df$to_col),
        cap  = edges_df$p_burn
    )
    source_arcs = data.frame(from = "S", to = ignition_nodes, cap = big_M)
    sink_arcs   = data.frame(from = target_nodes, to = "T",   cap = big_M)
    return_arc  = data.frame(from = "T", to = "S", cap = big_M)
    all_arcs = rbind(return_arc, source_arcs, sink_arcs, arcs)

    # Write NODES, ARCS+cap in AMPL .dat format
    f = file(outfile, open = "w")
    on.exit(close(f))
    writeLines("set NODES :=", f)
    writeLines(paste("  ", paste(all_nodes, collapse = " ")), f)
    writeLines(";", f)
    writeLines("param: ARCS: cap :=", f)
    for (k in seq_len(nrow(all_arcs)))
        writeLines(sprintf("  %-8s %-8s %.6f",
            all_arcs$from[k], all_arcs$to[k], all_arcs$cap[k]), f)
    writeLines(";", f)
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

**Model**

With a limited budget of B cells, we formulate a MILP to choose which cells to clear as firebreaks to minimise total weighted damage. Settlement cells have weight 10 and wetland cells weight 5. Ignition cells (row 21) and target cells cannot be cleared.

We try two propagation models that share the same budget constraint and objective:

$$\min \sum_{t \in \text{TARGETS}} w_t \cdot d_t \quad \text{s.t.} \quad \sum_{v} \text{clear}_v \leq B$$

**Model A (Binary propagation):** An edge $(u,v)$ is active if $p_{\text{burn}}(u,v) > \tau$. If a burning neighbour has an active edge to cell $v$ and $v$ is not cleared, then $v$ must burn:

$$\text{burn}_v \geq \text{burn}_u - \text{clear}_v \quad \forall (u,v) : p_{\text{burn}} > \tau$$

This is a worst-case model — every active edge carries fire with certainty. The threshold $\tau$ controls conservatism.

**Model B (Continuous reachability):** Instead of binary burn, we track a continuous reachability $x_v \in [0,1]$ for each cell. Fire propagates multiplicatively with compressed edge weights $w_{uv} = p_{\text{burn}}(u,v)^\alpha$:

$$x_v \geq w_{uv} \cdot x_u - \text{clear}_v \quad \forall (u,v) \in E$$

The compression $\alpha < 1$ prevents reachability from decaying below solver tolerance over many hops. The objective directly uses reachability at targets: $\min \sum_t w_t \cdot x_t$.

We solve both models at budget levels $B = 5, 8, 10, 13, 16, 20$, sweeping $\tau = 0, 0.05, \ldots, 0.20$ for Model A and $\alpha = 0.05, 0.1, 0.15, 0.2, 0.3$ for Model B.

**Code**

Model A — binary propagation MILP:

```ampl
set CELLS;
set IGNITION within CELLS;
set TARGETS within CELLS;
set EDGES within {CELLS, CELLS};

param weight {TARGETS};
param p_burn {EDGES};
param B >= 0, default 13;
param tau default 0;

var clear {CELLS} binary;
var burn {CELLS} binary;

minimize Total_Damage:
    sum {t in TARGETS} weight[t] * burn[t];

subject to Budget:
    sum {v in CELLS} clear[v] <= B;
subject to Ignition_Burns {u in IGNITION}:
    burn[u] = 1;
subject to Clear_Blocks {v in CELLS}:
    burn[v] <= 1 - clear[v];
subject to Propagation {(u,v) in EDGES: p_burn[u,v] > tau}:
    burn[v] >= burn[u] - clear[v];
```

Model B — continuous reachability MILP (key differences only):

```ampl
param alpha default 0.5;
param w {(u,v) in EDGES} = p_burn[u,v] ^ alpha;

var x {CELLS} >= 0, <= 1;

minimize Total_Damage:
    sum {t in TARGETS} weight[t] * x[t];

subject to Propagation {(u,v) in EDGES}:
    x[v] >= w[u,v] * x[u] - clear[v];
```

Both .run scripts use Gurobi with `mipgap=0.01 timelim=120`, then output objective, cleared cells, and damaged targets in CSV format. Model A outputs `burn[t]`; Model B outputs `x[t]`.

R function to generate the shared .dat file:

```r
#' Generate AMPL .dat file for Step 2 firebreak MILP
#'
#' @param edges_df  Data frame from spread_probabilities()
#' @param targets   Data frame with row, col, weight columns
#' @param outfile   Path to write the .dat file
#' @param landscape Integer matrix (21x21) of vegetation codes
#' @return Invisible NULL. Writes .dat file to disk.
generate_step2_dat = function(edges_df, targets, outfile, landscape) {
  node_name = function(r, c) paste0(r, "_", c)
  nr = nrow(landscape)
  nc = ncol(landscape)

  # --- 1. set CELLS: all vegetated (non-bare) cells ----
  cells = c()
  for (r in 1:nr) {
    for (c in 1:nc) {
      if (landscape[r, c] != 0) {
        cells = c(cells, node_name(r, c))
      }
    }
  }

  # --- 2. set IGNITION: vegetated cells in row 21 ----
  ignition = c()
  for (c in 1:nc) {
    if (landscape[nr, c] != 0) {
      ignition = c(ignition, node_name(nr, c))
    }
  }

  # --- 3. set TARGETS + param weight ----
  target_nodes  = node_name(targets$row, targets$col)
  target_weight = targets$weight

  # --- 4. set EDGES + param p_burn ----
  edge_from = node_name(edges_df$from_row, edges_df$from_col)
  edge_to   = node_name(edges_df$to_row,   edges_df$to_col)
  keep      = (edge_from %in% cells) & (edge_to %in% cells)
  edge_from = edge_from[keep]
  edge_to   = edge_to[keep]
  p_burn    = edges_df$p_burn[keep]

  # --- 5. Write .dat file ----
  con = file(outfile, open = "w")
  on.exit(close(con))
  wl = function(...) writeLines(paste0(...), con)

  wl("set CELLS :=")
  wl("  ", paste(cells, collapse = " "))
  wl(";")
  wl("")
  wl("set IGNITION :=")
  wl("  ", paste(ignition, collapse = " "))
  wl(";")
  wl("")
  wl("param: TARGETS: weight :=")
  for (i in seq_along(target_nodes)) {
    wl("  ", target_nodes[i], "  ", target_weight[i])
  }
  wl(";")
  wl("")
  wl("param: EDGES: p_burn :=")
  for (i in seq_along(edge_from)) {
    wl("  ", edge_from[i], "  ", edge_to[i], "  ",
       format(p_burn[i], digits = 8))
  }
  wl(";")
}
```

R function to call AMPL for one (model, B, param) combination:

```r
#' Run a single AMPL firebreak solve for either model
#'
#' @param model      Character "A" or "B"
#' @param B          Integer budget
#' @param param_name Character "tau" or "alpha"
#' @param param_val  Numeric parameter value
#' @return One-row data frame with damage, cleared cells, status
run_one = function(model, B, param_name, param_val) {
  tmp = tempfile(fileext = ".run")
  writeLines(c(
    "reset;",
    sprintf("model step2_model%s.mod;", model),
    "data step2.dat;",
    sprintf("let B := %d;", B),
    sprintf("let %s := %g;", param_name, param_val),
    sprintf("include step2_model%s.run;", model)
  ), tmp)
  output = system2(ampl_bin, tmp, stdout = TRUE, stderr = TRUE)
  unlink(tmp)
  # Parse SUMMARY and CLEARED lines from output
  parts = strsplit(grep("^SUMMARY,", output, value = TRUE), ",")[[1]]
  cleared = sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE))
  data.frame(model = model, B = B, param_name = param_name,
    param_value = param_val, damage = as.numeric(parts[2]),
    n_cleared = length(cleared),
    cleared_cells = paste(cleared, collapse = ";"),
    status = as.integer(parts[3]), stringsAsFactors = FALSE)
}
```

We sweep $B = 5, 8, 10, 13, 16, 20$ with $\tau = 0, 0.05, 0.10, 0.15, 0.20$ for Model A and $\alpha = 0.05, 0.1, 0.15, 0.2, 0.3$ for Model B. Marginal value is computed as $(\text{damage}(B) - \text{damage}(B+1)) / \Delta B$.

R function to validate plans with fire simulation:

```r
#' Run K fire simulations and return mean damage with 95% CI
#'
#' @param landscape  21x21 matrix (firebreaks applied as 0)
#' @param ignition   Two-column matrix (row, col)
#' @param targets    Data frame with row, col, weight
#' @param K          Number of replications
#' @return Named list: mean_damage, se, ci_lo, ci_hi
run_simulation = function(landscape, ignition, targets, K) {
  damages = numeric(K)
  for (k in seq_len(K)) {
    result = simulate_fire(landscape, ignition,
                           wind_speed = 0, wind_dir = 0)
    damages[k] = compute_damage(result$burned, targets)
  }
  mu = mean(damages)
  se = sd(damages) / sqrt(K)
  list(mean_damage = mu, se = se,
       ci_lo = mu - 1.96 * se, ci_hi = mu + 1.96 * se)
}
```

**Results**

We solved both models across all budget levels and validated the optimal plans with 1000 fire simulations (no wind, southern front ignition).

Model A shows a sharp threshold effect: at $\tau \leq 0.10$, fire deterministically reaches all targets (damage = 240); at $\tau \geq 0.105$, no fire reaches any target (damage = 0). The only interesting case is $\tau = 0$, $B = 13$, where the model places all 13 firebreaks to protect settlement (column 15 wall + row 7 seal), sacrificing wetland entirely (damage = 80).

Model B produces a smooth budget-damage curve. At $\alpha = 0.1$:

| Budget B | MILP damage | Sim. mean (K=1000) | 95% CI |
|----------|-------------|---------------------|--------|
| 0 (baseline) | — | 17.4 | [16.1, 18.7] |
| 5 | 38.2 | 10.7 | [8.5, 12.9] |
| 13 | 14.1 | 4.6 | [4.0, 5.1] |
| 16 | 12.0 | 2.7 | [1.9, 3.6] |
| 20 | 0 | 0 | [0, 0] |
| Corridor (17) | — | 1.6 | [1.1, 2.0] |

[INSERT: step2/step2_comparison.png]

**Figure 2:** Left: Model B damage vs budget for different $\alpha$ values. Right: MILP predictions vs simulator mean with 95% CI. Bottom: firebreak location maps for the current corridor, Model A ($B=16$, $\tau=0$), and Model B ($B=16$, $\alpha=0.1$).

**Discussion**

Model A is all-or-nothing: at $\tau = 0$ fire reaches every cell deterministically, so small budgets cannot help. Only $B \geq 13$ forms a complete barrier. Model B provides gradual tradeoffs — at $B = 13$ both models agree on the same plan (column 15 + row 7, protecting settlement), but at $B = 16$ Model B also protects wetland while Model A does not. Diminishing returns set in past $B = 13$, suggesting 13 cells is the minimum meaningful investment.

Both models overestimate damage — Model A predicts 80 but the simulator shows 4.6 for the same $B = 13$ plan. Yet the plan itself works well: the firebreak locations are in the right places. This shows a model can predict poorly but still produce a good plan. The parameter $\alpha$ requires experimentation: $\alpha = 0.5$ causes reachability to vanish over 15+ hops, while $\alpha = 0.1$ gives useful signal. The corridor (17 cells) achieves damage 1.6, slightly better than Model B at $B = 16$ (damage 2.7) but uses more cells.

These plans are optimised for no wind and southern ignition. Under different conditions the optimal locations could shift — west wind would push fire toward settlement faster. This is addressed in Step 3. If settlement and wetland had equal weights, the plan would split resources across both targets rather than concentrating on one side.

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
