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

Max-flow model (.mod):

```ampl
# Max-flow on fire spread network
# Step 1: Network bottleneck analysis
# Nodes: grid cells as "r_c" plus super-source "S" and "T"
# Edges: spread probabilities as capacities
# Return arc T -> S (unconstrained) for standard max-flow

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

AMPL script to solve max-flow and extract shadow prices (.run):

```ampl
option solver gurobi;

# --- Settlement ---
reset;
model step1_maxflow.mod;
data step1_settlement.dat;
solve;
printf "Settlement max flow = %g\n", TotalFlow;

printf "from,to,flow,capacity,shadow_price\n" > "shadow_prices_settlement.csv";
for {(i,j) in ARCS: i <> "T" or j <> "S"} {
    printf "%s,%s,%g,%g,%g\n",
        i, j, flow[i,j], cap[i,j], Capacity[i,j].dual
        >> "shadow_prices_settlement.csv";
}

# --- Wetland ---
reset;
model step1_maxflow.mod;
data step1_wetland.dat;
option solver gurobi;
solve;
printf "Wetland max flow = %g\n", TotalFlow;

printf "from,to,flow,capacity,shadow_price\n" > "shadow_prices_wetland.csv";
for {(i,j) in ARCS: i <> "T" or j <> "S"} {
    printf "%s,%s,%g,%g,%g\n",
        i, j, flow[i,j], cap[i,j], Capacity[i,j].dual
        >> "shadow_prices_wetland.csv";
}
```

R function to generate max-flow .dat file from spread probabilities:

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

Model A — binary propagation MILP (.mod):

```ampl
set CELLS;
set IGNITION within CELLS;
set TARGETS  within CELLS;
set EDGES    within {CELLS, CELLS};

param weight {TARGETS};
param p_burn {EDGES};
param B      >= 0, default 13;
param tau    default 0;

var clear {CELLS} binary;
var burn  {CELLS} binary;

minimize Total_Damage:
    sum {t in TARGETS} weight[t] * burn[t];

subject to Budget:
    sum {v in CELLS} clear[v] <= B;

subject to Ignition_Burns {u in IGNITION}:
    burn[u] = 1;

subject to Ignition_NoClear {u in IGNITION}:
    clear[u] = 0;

# Targets cannot be cleared: they must be protected by blocking fire upstream.
subject to Target_NoClear {t in TARGETS}:
    clear[t] = 0;

subject to Clear_Blocks {v in CELLS}:
    burn[v] <= 1 - clear[v];

# Only edges with p_burn > tau propagate fire (thresholded propagation).
subject to Propagation {(u, v) in EDGES: p_burn[u,v] > tau}:
    burn[v] >= burn[u] - clear[v];
```

Model B shares the same sets, budget, ignition, and target constraints as Model A. The differences are: `burn` is replaced by a continuous reachability `x` ∈ [0,1], edge weights are compressed as `w = p_burn^alpha` to prevent decay over many hops, and propagation is multiplicative rather than binary.

```ampl
param alpha default 0.5;
param w {(u,v) in EDGES} = p_burn[u,v] ^ alpha;

var x {CELLS} >= 0, <= 1;

minimize Total_Damage:
    sum {t in TARGETS} weight[t] * x[t];

subject to Ignition_Reach {u in IGNITION}:
    x[u] = 1;

subject to Clear_Blocks {v in CELLS}:
    x[v] <= 1 - clear[v];

subject to Propagation {(u, v) in EDGES}:
    x[v] >= w[u,v] * x[u] - clear[v];
```

Both models share a single `.run` skeleton:

```ampl
option solver gurobi;
option gurobi_options 'mipgap=0.01 timelim=120';
solve;

printf "SUMMARY,%f,%d\n", Total_Damage, solve_result_num;
for {v in CELLS: clear[v] > 0.5} printf "CLEARED,%s\n", v;
for {t in TARGETS: burn[t]  > 0.5} printf "DAMAGED,%s,%d,%f\n", t, weight[t], burn[t];
```

The Model B `.run` is identical except the damaged-target filter uses `x[t] > 0.01` instead of `burn[t] > 0.5`.

The shared `.dat` file is generated in R from the landscape matrix, the edge list returned by `spread_probabilities()` (Step 1), and `targets.csv`. The mapping is:

```r
#' Build AMPL .dat sets/params from Step-1 edge data and landscape matrix.
#'
#' @param edges_df  from spread_probabilities(): from_row, from_col, to_row, to_col, p_burn
#' @param targets   data frame: row, col, weight
#' @param landscape 21x21 integer matrix (0 = bare, non-zero = vegetated)
#' @return list with CELLS, IGNITION, TARGETS (+ weight), EDGES (+ p_burn);
#'   a thin writer then serialises each field as an AMPL set/param block.
generate_step2_dat = function(edges_df, targets, landscape) {
  node = function(r, c) paste0(r, "_", c)
  nr   = nrow(landscape)

  veg      = which(landscape != 0, arr.ind = TRUE)
  CELLS    = node(veg[, "row"], veg[, "col"])
  IGNITION = CELLS[veg[, "row"] == nr]                         # row 21

  TARGETS = data.frame(cell   = node(targets$row, targets$col),
                       weight = targets$weight)

  EDGES = data.frame(from   = node(edges_df$from_row, edges_df$from_col),
                     to     = node(edges_df$to_row,   edges_df$to_col),
                     p_burn = edges_df$p_burn)
  EDGES = EDGES[EDGES$from %in% CELLS & EDGES$to %in% CELLS, ]

  list(CELLS = CELLS, IGNITION = IGNITION, TARGETS = TARGETS, EDGES = EDGES)
}
```

Each field in that list is serialised to the shared `.dat` file as a single AMPL block, e.g. for `CELLS` and the `EDGES + p_burn` combined table:

```r
writeLines(c("set CELLS :=", paste(" ", CELLS, collapse = ""), ";"), con)
writeLines("param: EDGES: p_burn :=", con)
writeLines(sprintf("  %s  %s  %.8f", EDGES$from, EDGES$to, EDGES$p_burn), con)
writeLines(";", con)
```

`IGNITION` and `TARGETS + weight` follow the same pattern (`set` block and `param:` block respectively).

The AMPL solve is driven from R, writing a small driver `.run` per call so `B` and the swept parameter (`tau` or `alpha`) can be set externally:

```r
#' Solve one (model, B, parameter) combination and parse the CSV output.
#'
#' @param model      "A" or "B"
#' @param B          integer budget
#' @param param_name "tau" or "alpha"
#' @param param_val  numeric value for the swept parameter
#' @return one-row data frame: model, B, param_name, param_value,
#'   damage, n_cleared, cleared_cells, status
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

  parts   = strsplit(grep("^SUMMARY,", output, value = TRUE), ",")[[1]]
  cleared = sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE))

  data.frame(model = model, B = B,
             param_name = param_name, param_value = param_val,
             damage = as.numeric(parts[2]), n_cleared = length(cleared),
             cleared_cells = paste(cleared, collapse = ";"),
             status = as.integer(parts[3]))
}
```

We sweep $B = 5, 8, 10, 13, 16, 20$ with $\tau = 0, 0.05, 0.10, 0.15, 0.20$ for Model A and $\alpha = 0.05, 0.1, 0.15, 0.2, 0.3$ for Model B. Marginal value is reported as $(\text{damage}(B_j) - \text{damage}(B_{j+1})) / (B_{j+1} - B_j)$ for consecutive budgets inside each sweep.

Each optimal plan is validated by applying its cleared cells to the landscape and running $K = 1000$ replications of `simulate_fire()` (Step 3) under no-wind southern-front ignition. The baseline (no firebreaks) and current corridor (row 10, cols 3–19) are evaluated under the same protocol:

```r
#' Mean weighted damage + 95% CI for one firebreak plan.
#'
#' @param landscape 21x21 matrix with firebreaks applied (cleared cells -> 0)
#' @param ignition  two-column (row, col) ignition matrix
#' @param targets   data frame: row, col, weight
#' @param K         number of replications
validate_plan = function(landscape, ignition, targets, K = 1000) {
  set.seed(1)                                          # reproducible sampling
  damages = replicate(K, {
    sim = simulate_fire(landscape, ignition, wind_speed = 0, wind_dir = 0)
    compute_damage(sim$burned, targets)
  })
  mu = mean(damages); se = sd(damages) / sqrt(K)
  c(mean = mu, ci_lo = mu - 1.96 * se, ci_hi = mu + 1.96 * se)
}
```

**Results**

We solved both models across all budget levels and validated each optimal plan with K = 1000 fire simulations (no wind, southern-front ignition, `set.seed(1)`).

Model A exhibits a sharp threshold effect. At $\tau \leq 0.10$, every edge is active and fire deterministically reaches all targets (damage = 240); at $\tau \geq 0.105$, enough edges are filtered out that no fire reaches any target (damage = 0). The only useful case is $\tau = 0$ with $B = 13$, where the model places all 13 firebreaks along column 15 and row 7 to protect settlement, sacrificing wetland entirely (MILP damage = 80, simulated damage = 4.8).

Model B produces a smooth budget-damage curve. At $\alpha = 0.1$:

| Budget B | MILP damage | Sim. mean (K=1000) | 95% CI |
|----------|-------------|---------------------|--------|
| 0 (baseline) | — | 16.4 | [15.1, 17.7] |
| 5 | 38.2 | 10.1 | [9.1, 11.1] |
| 13 | 14.1 | 4.8 | [4.3, 5.4] |
| 16 | 12.0 | 3.0 | [2.5, 3.4] |
| 20 | 0 | 0 | [0, 0] |
| Corridor (17) | — | 1.6 | [1.1, 2.1] |

[INSERT: step2/step2_comparison.png]

**Figure 2:** Left: Model B damage vs budget for different $\alpha$ values. Centre: MILP predictions vs simulator mean with 95% CI. Right: firebreak location maps for the current corridor, Model A ($B=13$, $\tau=0$), and Model B ($B=13$, $\alpha=0.1$).

**Discussion**

Model A is all-or-nothing because at $\tau = 0$ every edge is active, so the propagation constraint forces every reachable cell to burn with certainty. Partial barriers are useless — fire routes around them through any remaining active edge. Only a complete wall ($B \geq 13$ along column 15 + row 7) blocks all paths to the settlement. Model B provides smooth tradeoffs: at $B = 13$ both models agree on the same plan (protecting settlement), but at $B = 16$ Model B also protects wetland while Model A still uses only 13 cells. The marginal value of additional cells drops sharply past $B = 13$, confirming 13 as the minimum meaningful investment.

Both models overestimate damage compared to the simulator — for instance, Model A predicts damage 80 at $B = 13$ but the simulator shows 4.8. The gap arises because the deterministic MILP treats every active edge as certain to propagate fire: if any burning neighbour has an active edge to cell $v$, then $v$ must burn. In practice, each edge only fires with probability $p_{\text{burn}}$, so most paths through the grid never fully activate in a single simulation run. The MILP therefore sees the worst-case superposition of all possible fire paths, while the simulator draws one stochastic realisation at a time. Despite this, the MILP plan is effective: the optimiser correctly identifies where fire pressure is highest, and the resulting firebreak locations block the critical corridors even though the predicted damage is far too high. This is a useful property — the model overpredicts damage but still produces a good plan, because the objective function ranks plans correctly even if the absolute values are wrong.

The parameter $\alpha$ requires experimentation: $\alpha = 0.5$ causes reachability to decay below solver tolerance over 15+ hops, while $\alpha = 0.1$ maintains useful signal across the full grid. The current corridor (17 cells) achieves mean damage 1.6, outperforming Model B at $B = 16$ (damage 3.0) but using one more cell.

These plans are optimised for no-wind, southern-front ignition. Under wind, the optimal firebreak locations would shift — for example, westerly wind would push fire toward the settlement faster, potentially requiring a wider barrier on that side. Step 3 evaluates plan robustness under uncertain wind conditions.

---

### Step 3: Evaluation Under Uncertain Wind

**Model**

We evaluate the proposed plan ($B = 13$) under stochastic wind and random ignition. Wind speed $V \sim \text{Weibull}(k=2, A=7)$, sampled via inverse transform. Wind direction follows a mixture of two von Mises:

$$f_\Theta(\theta) = 0.7 \cdot \frac{e^{3\cos(\theta - 1.5\pi)}}{2\pi\, I_0(3)} + 0.3 \cdot \frac{e^{2\cos(\theta - 0.75\pi)}}{2\pi\, I_0(2)}$$

The dominant mode ($\mu_1 = 1.5\pi$, westerly) pushes fire toward the settlement; the secondary mode ($\mu_2 = 0.75\pi$, south-easterly) toward the wetland. The von Mises CDF has no closed form, so we use acceptance-rejection with Uniform$(0, 2\pi)$ proposal, envelope $M = \max_\theta f_\Theta(\theta)$, and acceptance condition $U < f_\Theta(\theta)/M$.

Each of $N$ scenarios samples $(V, \Theta, \text{ignition})$, runs the simulator, and records damage. The three plans are evaluated on the same sampled scenarios, so any difference in risk measures reflects the firebreaks, not the draws. We report $\mathbb{E}[D]$ with a 95% $t$-interval and $\text{CVaR}_{0.9} = \mathbb{E}[D \mid D \geq q_{0.9}]$ with a bootstrap CI.

**Code**

Three files implement the evaluation: `step3_samplers.R`, `step3_monte_carlo.R`, and `step3_run.R`. `simulate_fire()`, `set_firebreaks()`, and `compute_damage()` are from the provided `fire-simulator.r` and used unchanged.

Samplers (`step3_samplers.R`):

```r
#' Implements the wind-direction density f_Theta(theta) defined in the Model.
dwind_dir = function(theta) {
  p1  = 0.7;   mu1 = 1.5 * pi;  kappa1 = 3
  p2  = 0.3;   mu2 = 0.75 * pi; kappa2 = 2

  vm1 = exp(kappa1 * cos(theta - mu1)) / (2 * pi * besselI(kappa1, nu = 0))
  vm2 = exp(kappa2 * cos(theta - mu2)) / (2 * pi * besselI(kappa2, nu = 0))

  p1 * vm1 + p2 * vm2
}

#' Samples wind direction by acceptance-rejection with Uniform(0, 2*pi) proposal.
#' Envelope M = max f_Theta(theta) obtained by optimise().
#' Returns samples plus the empirical and theoretical acceptance rates.
sample_wind_dir = function(n) {
  opt = optimise(dwind_dir, interval = c(0, 2 * pi), maximum = TRUE)
  M   = opt$objective

  samples      = numeric(n)
  n_accepted   = 0L
  n_total      = 0L
  total_accept = 0L

  while (n_accepted < n) {
    batch_size = max(n - n_accepted, 256L)
    theta = runif(batch_size, min = 0, max = 2 * pi)
    u     = runif(batch_size)

    accept       = u < dwind_dir(theta) / M
    n_total      = n_total + batch_size
    n_new        = sum(accept)
    total_accept = total_accept + n_new
    if (n_new == 0L) next

    n_take = min(n_new, n - n_accepted)
    samples[(n_accepted + 1):(n_accepted + n_take)] = theta[accept][1:n_take]
    n_accepted = n_accepted + n_take
  }

  list(samples = samples, acceptance_rate = total_accept / n_total,
       theoretical_rate = 1 / (2 * pi * M), M = M)
}

#' Samples wind speed (km/h) from Weibull(shape = 2, scale = 7).
sample_wind_speed = function(n) rweibull(n, shape = 2, scale = 7)

#' Samples ignition cells uniformly over the vegetated cells of `landscape`.
sample_ignition = function(landscape, n) {
  veg_idx = which(landscape != 0, arr.ind = TRUE)
  chosen  = veg_idx[sample(nrow(veg_idx), size = n, replace = TRUE), , drop = FALSE]
  colnames(chosen) = c("row", "col")
  chosen
}
```

Monte Carlo evaluation and risk measures (`step3_monte_carlo.R`):

```r
#' Evaluates all plans on a pre-sampled batch of (wind_dir, wind_speed, ignition)
#' scenarios, supplied by the caller so the same draws feed both diagnostics
#' and the MC. Ignition on a cleared cell yields zero damage (no fire starts).
run_mc_evaluation = function(landscape, targets, plans,
                             wind_dir, wind_speed, ignition) {
  N = length(wind_dir)
  stopifnot(length(wind_speed) == N, nrow(ignition) == N)
  results = list()
  for (plan_name in names(plans)) {
    ls_plan = landscape
    breaks  = plans[[plan_name]]
    if (!is.null(breaks)) ls_plan = set_firebreaks(ls_plan, breaks)

    damages = numeric(N)
    for (i in seq_len(N)) {
      ir = ignition[i, "row"]; ic = ignition[i, "col"]
      if (ls_plan[ir, ic] == 0) {
        damages[i] = 0
      } else {
        result     = simulate_fire(ls_plan, ignition[i, , drop = FALSE],
                                   wind_speed[i], wind_dir[i])
        damages[i] = compute_damage(result$burned, targets)
      }
    }
    results[[plan_name]] = data.frame(scenario = seq_len(N), damage = damages,
                                      wind_speed = wind_speed, wind_dir = wind_dir,
                                      ign_row = ignition[, "row"],
                                      ign_col = ignition[, "col"])
  }
  results
}

#' Computes mean damage, a 95% t-interval, and bootstrap CI for CVaR_alpha.
compute_risk_measures = function(damages, alpha = 0.9, n_boot = 2000) {
  n     = length(damages)
  mu    = mean(damages)
  se_mu = sd(damages) / sqrt(n)
  ci_mu = mu + c(-1, 1) * qt(0.975, df = n - 1) * se_mu

  cvar_point = function(x) {
    threshold = quantile(x, probs = alpha)
    mean(x[x >= threshold])
  }
  cvar_est  = cvar_point(damages)
  boot_cvar = replicate(n_boot,
                        cvar_point(damages[sample.int(n, replace = TRUE)]))
  ci_cvar   = quantile(boot_cvar, probs = c(0.025, 0.975))

  data.frame(measure  = c("mean_damage", paste0("CVaR_", alpha)),
             estimate = c(mu, cvar_est),
             ci_lower = c(ci_mu[1], ci_cvar[[1]]),
             ci_upper = c(ci_mu[2], ci_cvar[[2]]))
}
```

Workflow (`step3_run.R`). Using the landscape matrix and target table defined earlier from `landscape.csv` and `targets.csv`:

```r
set.seed(1)
N = 5000

## Proposed plan: Step 2 Model B, B=13, alpha=0.1 (L-shape at col 15 / row 7).
proposed = rbind(cbind(row = 1:7,    col = 15L),
                 cbind(row = 7L,     col = 16:21))
corridor = cbind(row = rep(10L, 17), col = 3:19)
plans    = list(baseline = NULL, corridor = corridor, proposed = proposed)

## Pre-sample N scenarios ONCE on the set.seed(1) stream.
dir_draws  = sample_wind_dir(N)
wind_dir   = dir_draws$samples
wind_speed = sample_wind_speed(N)
ignition   = sample_ignition(landscape, N)

## Diagnostics reuse the SAME draws — no extra RNG consumed.
pwind_dir = function(x) {
  sapply(x, function(t) integrate(dwind_dir, 0, t)$value)
}
ks_p = ks.test(pwind_dir(wind_dir), "punif")$p.value

cat(sprintf("M = %.3f   acceptance rate = %.3f   theoretical = %.3f   KS p = %.2f\n",
            dir_draws$M, dir_draws$acceptance_rate,
            dir_draws$theoretical_rate, ks_p))

## Monte Carlo evaluation on those same scenarios.
results    = run_mc_evaluation(landscape, targets, plans,
                               wind_dir, wind_speed, ignition)
risk_table = do.call(rbind, lapply(names(results), function(nm) {
  cbind(plan = nm, compute_risk_measures(results[[nm]]$damage))
}))
```

**Results**

The sampler achieves $M = 0.464$, acceptance rate 33.8% (theoretical $1/(2\pi M) = 34.3\%$). KS test $p = 0.73$ (Figure 3).

[INSERT: figures/fig1_sampler_validation.png]

**Figure 3:** Wind direction samples vs mixture density (left); wind speed samples vs Weibull(2, 7) (right).

$N = 5{,}000$ scenarios (`set.seed(1)`):

| Plan | Cells | $\mathbb{E}[D]$ | 95% CI | $\text{CVaR}_{0.9}$ | 95% CI |
|------|-------|---------|--------|----------|--------|
| No firebreaks | 0 | 13.77 | [13.19, 14.35] | 58.94 | [57.48, 66.21] |
| Corridor | 17 | 4.09 | [3.77, 4.42] | 23.83 | [22.61, 34.28] |
| Proposed | 13 | 4.34 | [4.08, 4.61] | 26.38 | [25.31, 27.50] |

[INSERT: figures/fig2_damage_densities.png]

**Figure 4:** Damage distributions. Dashed = mean, solid = $\text{CVaR}_{0.9}$.

[INSERT: figures/fig3_proposed_diagnostics.png]

**Figure 5:** Proposed plan: damage vs wind speed by compass direction (left); mean damage by ignition cell (right).

$N = 1{,}000$ pilot gave $\mathbb{E}[D] = 4.69$ vs 4.34, confirming stability.

**Discussion**

Both plans cut $\mathbb{E}[D]$ from 13.8 to ~4, but $\text{CVaR}_{0.9} \approx 25$ (6$\times$ the mean) reveals heavy tail risk. Corridor and proposed CIs overlap, so the difference is not significant — but proposed uses 4 fewer cells, giving better per-cell efficiency.

The worst 10% are driven by ignition near targets (bypassing firebreaks) and westerly wind pushing fire toward the settlement (Figure 5). The plan handles westerly wind well since its firebreaks block settlement-bound paths, but leaves the wetland exposed to south-easterly wind.

CVaR is a more informative criterion than $\mathbb{E}[D]$: the mean is diluted by many zero-damage scenarios where fire never reaches targets. A CVaR-aware design would place firebreaks closer to targets to guard against nearby ignitions.

Step 2's MILP predicted damage 14.1; the stochastic mean is 4.34. The plan remains effective because the MILP identified the structurally important corridors, which wind does not fundamentally alter.

---

### Step 4: Reactive Response

**Model**

Step 4 splits the combined budget $B + R = 17$ into two stages. The $B = 13$ permanent cells from Step 3 are fixed before the fire season. After wind conditions $(V, \Theta)$ are revealed — but before ignition — up to $R = 4$ additional cells may be cleared reactively.

The reactive cells are chosen by a per-scenario MILP that extends Model B. The permanent cells enter as fixed data (set PERM), and the only decision variables are binary $z_v$ for each candidate cell:

$$\min \sum_{t \in \text{TARGETS}} w_t \, x_t \quad \text{s.t.} \quad \sum_v z_v \le R$$

with $w_{uv} = p_{\text{burn}}(u \to v;\, V, \Theta)^\alpha$ computed from the revealed wind using `spread_probabilities()`. The propagation and blocking constraints from Model B apply with $\text{clear}[v]$ replaced by $\text{perm}[v] + z[v]$; edges involving permanent cells are filtered out. The compression exponent $\alpha = 0.1$ is unchanged from Step 2.

The true ignition location is unknown at both stages. Following the brief, the reactive MILP uses a southern-front proxy: all row-21 cells have $x_u = 1$, consistent with Step 2. This is conservative; Q6 discusses its limitations.

We compare four plans under the same $N = 1000$ scenarios sampled with `set.seed(1)`: (1) no firebreaks, (2) $B = 13$ permanent only, (3) $B + R = 17$ all-permanent, and (4) $B = 13$ permanent + $R = 4$ reactive. Within each scenario, `set.seed(1000 + i)` aligns the simulator's starting RNG state across plans, reducing noise in per-scenario comparisons without creating exact pathwise coupling (different firebreak layouts consume `runif` at different rates).

**Code**

Step 4 extends Model B with a reactive MILP, a wind-aware `.dat` builder (`edges_wind.R`, `build_step4_dat.R`), and a shared-scenario batch evaluator. Samplers and risk measures from Step 3 are reused.

Reactive MILP — Step 4 constraints beyond Model B (`step4_reactive.mod`):

```
set PERM within CELLS;
param R >= 0;
var react {CELLS} binary;

check {u in IGNITION}: u not in PERM;
check {t in TARGETS}:  t not in PERM;

subject to Reactive_Budget:
    sum {v in CELLS} react[v] <= R;
subject to Ignition_NoReact {u in IGNITION}: react[u] = 0;
subject to Target_NoReact {t in TARGETS}:    react[t] = 0;
subject to Perm_Reach {v in PERM}:           x[v] = 0;
subject to Perm_NoReact {v in PERM}:         react[v] = 0;
subject to React_Blocks {v in CELLS diff PERM}:
    x[v] <= 1 - react[v];
subject to Propagation {(u,v) in EDGES: u not in PERM and v not in PERM}:
    x[v] >= w[u,v] * x[u] - react[v];
```

The `.dat` builder (`build_step4_dat.R`) recomputes edge weights for the revealed wind via `spread_probabilities()` and writes PERM, $w$, and $R$ alongside the standard sets:

```r
wl("set PERM :=")
wl("  ", paste(perm_nodes, collapse = " "))
wl(";")

wl("param: EDGES: w :=")
for (i in seq_along(edge_from))
  wl("  ", edge_from[i], "  ", edge_to[i], "  ", format(w_vals[i], digits = 8))
wl(";")

wl("param R := ", R, ";")
```

Evaluation loop — core of `evaluate_all_plans.R`, showing the shared-scenario paired comparison with per-scenario RNG alignment:

```r
for (i in seq_len(N)) {
  sim_seed = 1000L + i

  set.seed(sim_seed)
  res1 = simulate_fire(land_none, ign, V, dir)
  d1   = compute_damage(res1$burned, targets)       # plan 1: none

  set.seed(sim_seed)
  res2 = simulate_fire(land_13, ign, V, dir)
  d2   = compute_damage(res2$burned, targets)       # plan 2: perm_13

  set.seed(sim_seed)
  res3 = simulate_fire(land_17, ign, V, dir)
  d3   = compute_damage(res3$burned, targets)       # plan 3: perm_17

  z    = solve_reactive(perm_13, V, dir, alpha, R)  # MILP: wind only
  set.seed(sim_seed)
  res4 = simulate_fire(set_firebreaks(land_13, z), ign, V, dir)
  d4   = compute_damage(res4$burned, targets)       # plan 4: reactive
}
```

Each scenario resets the simulator's RNG to `1000 + i` before each plan, aligning starting states without exact pathwise coupling. `solve_reactive` is conditioned on $(V, \Theta)$ and the permanent plan only; the realised ignition enters only in the subsequent simulator evaluation. For sensitivity, we reuse the same evaluator with alternative permanent plans and reactive budgets.

**Results**

We evaluate four plans on the same $N = 1{,}000$ scenarios (Table 1). $N = 1{,}000$ balances precision against runtime: each scenario requires a reactive MILP solve, so Step 4's per-scenario cost is roughly two orders of magnitude larger than Step 3's simulator-only evaluation. At this $N$, point estimates are stably ordered, and the no-firebreak and $B = 13$ plans are clearly worse than the two 17-cell plans on both mean and tail risk. However, the two 17-cell plans overlap on both mean and CVaR intervals, and $B = 13$ still overlaps them on CVaR, so we do not claim a statistically resolved ordering from intervals alone.

| Plan | $\mathbb{E}[D]$ (95% $t$-CI) | $\text{CVaR}_{0.9}$ (95% bootstrap CI) |
|---|---|---|
| No firebreaks | $14.27 \pm 1.32$ | $64.47\ [55.10,\ 69.25]$ |
| $B = 13$ permanent | $4.84 \pm 0.60$ | $28.45\ [24.57,\ 30.55]$ |
| $B + R = 17$ all-permanent | $3.54 \pm 0.54$ | $25.95\ [18.42,\ 28.37]$ |
| $B = 13$ + reactive ($R = 4$) | $3.34 \pm 0.53$ | $19.90\ [17.75,\ 27.10]$ |

**Table 1:** Risk measures for the four plans over $N = 1{,}000$ shared scenarios. Mean CI uses the $t$-approximation; CVaR CI is bootstrap with $2{,}000$ resamples.

The scenario-wise gap $D_{\text{perm\_13}} - D_{\text{reactive}}$ has mean $+1.49$; reactive strictly reduces damage in $16.8\%$ of scenarios and strictly increases it in $6.9\%$ (plans (2) and (4) share scenario-level RNG starting states but not exact pathwise realisations, so some drift is expected; Figure 6a). Gurobi reached the $1\%$ target MIP gap on $960$ of $969$ reactive solves; median final gap was $0.97\%$. The remaining $31/1{,}000$ scenarios had ignition on a permanent firebreak, trivially giving zero damage.

*Where reactive effort matters most.* Stratifying the gap by covariates (Figure 6b) gives a partial answer. Ignition location shows the clearest signal: for ignitions in rows 1–7 (near the targets) the mean gap is $0.37$ and $P(\text{gap} > 0) = 3.2\%$; in rows 8–14 the mean gap rises to $2.68$ ($P(\text{gap} > 0) = 22.3\%$); in rows 15–21 (near the MILP's southern-front proxy) the mean gap is $1.38$. Wind-speed terciles do not show a monotone trend (mean gaps $1.87 / 0.92 / 1.70$ from low to high), so the evidence does not support the claim that reactive matters most at extreme wind speeds. Against `perm_17` specifically, the mean damage difference is essentially zero under westerly wind ($\Delta = 0.03$, $n = 714$) and larger under easterly wind ($\Delta = 0.65$, $n = 286$). These stratifications primarily reflect where the southern-front proxy aligns or misaligns with the realised fire geometry; they do not yet isolate a pure "marginal value of wind information" signal.

**Figure 6:** (a) Distribution of the scenario-wise reactive effort gap $D_{\text{perm\_13}} - D_{\text{reactive}}$ with mean (solid) and zero (dashed). (b) Mean gap split by ignition row: reactive improves most when ignition lies in the middle rows (8–14), less at either end.

*Reactive cell frequency* (Figure 7). Four cells dominate: $(11, 4)$, $(11, 5)$, $(11, 6)$ are each selected in $> 96\%$ of scenarios and $(15, 5)$ in $73.7\%$ — exactly the four cells the no-wind MILP at $B = 17$ adds on top of $B = 13$. Of the $3{,}876$ total reactive selections, $99.6\%$ are closer (Chebyshev distance) to the wetland than to the settlement; $0\%$ are closer to the settlement.

**Figure 7:** Spatial frequency of reactive cell selection across $N = 1{,}000$ scenarios. Colour intensity = number of scenarios selecting that cell. Black outlines mark the $B = 13$ permanent plan; red and blue outlines mark the settlement and wetland targets. Reactive effort concentrates on the four cells the no-wind $B = 17$ MILP adds to $B = 13$.

*$R$ sensitivity* (Figure 8). At $R = 0$ plan (4) reduces to `perm_13` (E $= 3.70$, $N = 200$); $R = 2$ lowers the mean to $3.47$ but is anomalous in $\text{CVaR}_{0.9}$ ($24.05$) relative to the broader trend; $R = 4$ gives E $= 3.34$, CVaR $= 19.90$ ($N = 1{,}000$). The observed $R = 0 \to 4$ improvement is clear. At larger $R$, point estimates continue to fall: $R = 5$ gives E $= 2.57$, CVaR $= 17.69$ ($131/193$ solves proved optimal, median MIP gap $1\%$); $R = 6$ gives E $= 1.25$, CVaR $= 10.87$; $R = 8$ gives E $= 1.47$, CVaR $= 10.93$. Solver quality degrades as $R$ grows — at $R = 6$ only $5/193$ solves were proved optimal (median MIP gap $27\%$), so that point alone would be suggestive at best; the intermediate $R = 5$ point is consistent with the downward trend, which makes the drop from $R = 4$ to $R \ge 5$ less likely to be explained entirely by solver noise. The $R = 6 \to R = 8$ non-monotonicity is within $N = 200$ sampling noise and the two points should be treated as indistinguishable. Overall the evidence is consistent with diminishing returns around $R = 6$–$8$, but the exact elbow is not pinned down by the current runs.

**Figure 8:** $R$ sensitivity of plan (4). (a) Mean damage and (b) $\text{CVaR}_{0.9}$ as functions of the reactive budget $R \in \{0, 2, 4, 5, 6, 8\}$, with `perm_13` and `perm_17` as horizontal references. The $R = 6$ point is circled because its MIP gap is large and the result is not certified optimal.

**Discussion**

**Q1.** A fully optimal approach is a two-stage stochastic program with non-anticipative recourse, where the reactive policy $z(V, \Theta)$ depends on the revealed wind but not on the unobserved ignition $I$:
$$\min_{\text{perm},\ z(\cdot)}\ \mathbb{E}_{V,\Theta}\bigl[\ \mathbb{E}_{I}\bigl[\ D(\text{perm},\ z(V,\Theta);\ V, \Theta, I)\ \bigr]\ \bigr].$$
This requires scenario-based decomposition over the joint $(V, \Theta, I)$ distribution, which is computationally prohibitive on a $21 \times 21$ grid at a feasible budget. Our Step 4 surrogate replaces the inner expectation with a southern-front proxy and solves one MILP per sampled $(V, \Theta)$.

**Q2.** Access to realised wind weakly enlarges the reactive policy class, so plan (4) should match or beat plan (3) in expectation; the surrogate does not guarantee exact sample-level dominance. The observed means are consistent ($3.34$ vs $3.54$; $\text{CVaR}_{0.9}$ $19.90$ vs $25.95$), driven mostly by the minority easterly-wind regime ($\Delta = 0.65$, $n = 286$) — under dominant westerly wind the two plans are indistinguishable on average ($\Delta = 0.03$, $n = 714$). The $2.4\%$ of scenarios where reactive exceeds `perm_17` are consistent with the ignition proxy (Q6) and residual simulator-RNG drift rather than a systematic modelling failure.

**Q3.** Gurobi reached the $1\%$ target gap on $99.1\%$ of reactive solves, so the optimality gap does not materially affect the reactive cell choices in our sample. Separately, the $B = 17$ all-permanent MILP reported a $25\%$ gap at the $120$ s limit, but a re-solve at $1{,}800$ s confirmed the same incumbent is globally optimal. A large reported gap here reflects a slow LP lower bound rather than a suboptimal incumbent.

**Q4.** The reactive MILP concentrates its four cells on $(11,4), (11,5), (11,6), (15,5)$ in most scenarios, matching the no-wind MILP's $B = 17$ additions to $B = 13$. Wind induces only small secondary shifts. This suggests the permanent plan leaves one specific structural gap on the wetland side that is almost always worth closing, rather than the reactive stage exploiting fine-grained wind-specific adjustments.

**Q5.** The settlement weight (10) is twice the wetland weight (5), yet $99.6\%$ of reactive-cell selections are closer to the wetland. This is consistent with the permanent plan already largely saturating settlement protection — the $B = 13$ L-shape blocks the main northward corridor — so the marginal damage reduction from additional settlement-side cells is small and reactive redirects $R = 4$ to the wetland. Under this reading, weight asymmetry is absorbed by the permanent plan rather than expressed in the reactive layer.

**Q6.** The reactive MILP is conditioned on $(V, \Theta)$ only; ignition $I$ is unobserved, so the MILP requires a surrogate. We follow Step 2 and use a southern-front proxy ($x_u = 1$ for row 21), one approximation among several. The stratified gaps in Results are consistent with the expected degradation as ignition moves away from the proxy front: for ignitions in the northern third, reactive rarely improves over `perm_13` ($P(\text{gap} > 0) = 3.2\%$); in the middle rows the mean gap is largest ($2.68$). A deterministic alternative (a wind grid weighted by the mixture density) would improve reproducibility but lose the bootstrap-based CVaR uncertainty quantification, so we retain Monte Carlo.

**Q7.** Re-solving Step 2 Model B at $B = 12$ gives a 12-cell permanent plan (`perm_12`) equal to `perm_13` with cell $(1, 15)$ removed. Evaluated on $N = 200$ scenarios, the $(B = 12, R = 5)$ split gives reactive mean $3.15$ and $\text{CVaR}_{0.9} = 18.45$, compared with $3.34$ and $19.90$ for the main $(B = 13, R = 4)$ split. The two splits are within each other's CI, suggesting the allocation is not sensitive within this range: trading one permanent cell for one reactive cell does not meaningfully change risk outcomes under our surrogate.

---

## Discussion

(to write after all steps — synthesise across steps)

---

## AI Declaration

(to write last)
