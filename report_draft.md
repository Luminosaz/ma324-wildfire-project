# Submission

## Current Solution

The  current plan clears 17 cells along row 10 (column 3-19), forming a single east-west barrier between the southern ignition front and the northern targets.

Rationale that might support it. A mid-grid horizontal barrier is a simple static plan: it is installed once, does not rely on weekly forecast information, and treats both targets uniformly without requiring an ex ante decision about how protection should be prioritised between them Step 3's stochastic evaluation ($N = 5000$）also shows that the corridor is not ineffective: in the same regime it reduces mean damage from $\mathbb{E}[D] = 13.77$ with no firebreak to $\mathbb{E}[D] = 4.09$. The issue is therefore not that the corridor fails completely, but that later analysis identified better-supported alternatives.



Problem identifed.





Misalignment with Step 1 bottlenecks. Shadow prices (Figure 1) peak at rows 6-7, adjacent to each target; row 10 is well south of these peaks, so the corridor intercepts fire far from where Step 1's flow network pinches.



Insensitive to the 10:5 weight asymmetry. A symmetric full-width barrier spreads the budget uniformly across both targets and does not concentrate protection on the settlement, whose damage weighs twice as much in the objective.



Static under a weekly forecast regime. The corridor is fixed before the season, so the same layout is deployed regardless of the realised wind.



In addition, the corridor leaves columns 1-2 and 20-21 uncleared; we note this as a structural gap, not a quantified loss.



Stochastic performance ( Step 3, $N = 5000$). The same evaluation also exposes substantial tail risk: ${CVaR}_{0.9} = 23.83$ — tail-to mean ratio$\approx 5.8$. Worst scenarios combine target-adjacent ignitions with strong wind.



## Proposal

We recommend replacing the current 17-cell row-10 corridor with a two-stage firebreak plan at the same 17-cell footprint cap.



Before the season, install 13 permanent firebreak cells in an L-shape: column 15, rows 1-7 , and row 7, columns 16-21. After wind speed and direction are revealed, but before ignition, choose up to four additional reactive cells using the Step 4 recourse rule. The reactive cells are not fixed in advance; they are selected for the realised wind scenario, subject to the same legality constraints as the permanent plan.



This proposal addresses the main weaknesses of the current corridor. First, the permanent L-shape spans row 1-7 and therefore covers the row 6-7 band where Step 1 shadow prices are strongest, rather than at row 10 where the corridor sits well south of the bottleneck. Second, it anchors the barrier to the northern boundary ( at column 15, row 1) and the eastern boundary (at row 7, column 21), closing the approaches that a mid-grid corridor leaves open on its flanks. Third, it places the permanent budget near the higher-value settlement while retaining a reactive layer that can respond to revealed wind, instead of committing all 17 cells before season.



This two-stage plan is best read as the preferred evaluated extension of the certified $B = 13$ L-shape scaffold, rather than as a globally resolved best 17-cell policy across all candidate permanent plans. The current corridor is not ineffecive, but the later analysis shows that the row-10 geometry is structurally misaligned with the main Step 1 bottlenecks and that substantial tail risk remains under stochastic wind and ignition. The fuller quantitative justification, including the comparative Monte Carlo evaluation and the limits of the corridor-versus-two-stage comparison, is given later in the report.



The proposal is feasible. It respects the overall 17-cell footprint cap, using 13 permanent cells plus at most four reactive cells. All permanent cells are vegetated and legally clearable, and ignition- row and target cells are excludes from both stages by construction. The reactive stage is conditional on revealed wind only, not on the realised ignition cell, so the information structure is consistent with the Step 4 specification.





# Scaffolding

## Step 1: Network Bottleneck Model

### Model

For Step 1, we model the fire grid as a network flow problem. Each vegetated cell on 21x21 grid becomes a node, and we draw directed edges between adjacent cells (including diagonals, so 8 directions). The capacity of each edge is set to the spread probability p_brn from the simulator with no wind.



We add a source node S that connects to all cells in row 21 (the south edge where fire starts), and a sink node T that connects from the target cells. There is also a return arc from T back to S so we can use the standard max-flow formulation.



The LP maximises total flow from S to T:



$$\begin{aligned} \max \quad & \text{flow}(T, S) \\ \text{s.t.} \quad & \sum_{(i,v) \in E} \text{flow}(i,v) = \sum_{(v,j) \in E} \text{flow}(v,j) \quad \forall \, v \in V \\ & \text{flow}(i,j) \leq \text{cap}(i,j) \quad \forall \, (i,j) \in E \end{aligned}
$$


We write the capacity constraint as a named constraint in AMPL so we can read shadow prices with .dual. We run this twice - once for settlement as the target, once for wetland.



To figure out which cells are most important to clear, we look at the shadow prices. Since clearing a cell removes all edges going into it, we sum up the shadow prices of all incoming edges for each cell. A higher sum means that cell is a bigger bottleneck.



### Code

Max-flow model(step1_maxflow.mod):

```
set NODES;
set ARCS within {NODES, NODES};
param cap {ARCS} >= 0;
var flow {(i,j) in ARCS} >= 0;

# Return arc T -> S closes the circulation so max-flow = flow["T","S"]
maximize TotalFlow: flow["T","S"];

subject to FlowBalance {v in NODES}:
    sum {(i,v) in ARCS} flow[i,v] - sum {(v,j) in ARCS} flow[v,j] = 0;

# Named (not bound) so Capacity[i,j].dual gives shadow prices
subject to Capacity {(i,j) in ARCS: i <> "T" or j <> "S"}:
    flow[i,j] <= cap[i,j];
```

Build the .dat file for both targets under no wind:

```
source("fire-simulator.r")

landscape = as.matrix(read.csv("landscape.csv", header = FALSE))
targets   = read.csv("targets.csv")
edges     = spread_probabilities(landscape, wind_speed = 0, wind_dir = 0)

settlement_cells = targets[targets$weight == 10, c("row", "col")]
wetland_cells    = targets[targets$weight == 5,  c("row", "col")]

#' Generate AMPL .dat file for the max-flow model
#'
#' Builds the network (super-source -> ignition row, spread edges,
#' targets -> super-sink, return arc) in AMPL-readable format
#' matching step1_maxflow.mod.
#'
#' @param edges_df data frame from spread_probabilities() with
#'   from_row, from_col, to_row, to_col, p_burn
#' @param target_cells data frame with row, col columns
#' @param outfile path for the output .dat file
#' @param ignition_row fire source row (default 21 = southern front)
#' @param big_M large capacity for source/sink/return arcs
#' @return Writes an AMPL .dat file to `outfile`
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
    if (length(target_nodes)   == 0) stop("No target nodes found")

    arcs = data.frame(
        from = node_name(edges_df$from_row, edges_df$from_col),
        to   = node_name(edges_df$to_row,   edges_df$to_col),
        cap  = edges_df$p_burn
    )

    all_arcs = rbind(
        data.frame(from = "T", to = "S", cap = big_M),
        data.frame(from = "S", to = ignition_nodes, cap = big_M),
        data.frame(from = target_nodes, to = "T",   cap = big_M),
        arcs
    )

    f = file(outfile, open = "w")
    on.exit(close(f))

    writeLines("set NODES :=", f)
    for (chunk in split(all_nodes, ceiling(seq_along(all_nodes) / 10))) {
        writeLines(paste("   ", paste(chunk, collapse = "  ")), f)
    }
    writeLines(";", f)
    writeLines("", f)

    writeLines("param: ARCS: cap :=", f)
    for (k in seq_len(nrow(all_arcs))) {
        writeLines(sprintf("    %-8s %-8s  %.6f",
                           all_arcs$from[k], all_arcs$to[k], all_arcs$cap[k]), f)
    }
    writeLines(";", f)
}

generate_maxflow_dat(edges, settlement_cells, "step1_settlement.dat")
generate_maxflow_dat(edges, wetland_cells,    "step1_wetland.dat")
```



Solve script (step1_maxflow.run) — loops over both targets and writes one CSV per scenario:

```
option solver gurobi;

for {scen in {"settlement", "wetland"}} {
    reset;
    model step1_maxflow.mod;
    data ("step1_" & scen & ".dat");
    solve;

    printf "SUMMARY,%s,%g,%s,%d\n",
        scen, TotalFlow, solve_result, solve_result_num;

    printf "from,to,flow,capacity,shadow_price\n"
        > ("shadow_prices_" & scen & ".csv");

    for {(i,j) in ARCS: i <> "T" or j <> "S"} {
        printf "%s,%s,%g,%g,%g\n",
            i, j, flow[i,j], cap[i,j], Capacity[i,j].dual
            >> ("shadow_prices_" & scen & ".csv");
    }
}
```



Aggregate per-cell bottleneck score from the shadow-price CSVs:

```
#' Aggregate incoming-arc shadow prices to a cell-level bottleneck score
#'
#' @param csv_path path to shadow_prices_*.csv written by AMPL
#' @return data frame with row, col, and aggregated incoming dual
aggregate_bottleneck = function(csv_path) {
    sp = read.csv(csv_path)

    # Keep only real cell destinations; drop artificial source/sink arcs.
    sp = sp[grepl("^[0-9]+_[0-9]+$", sp$to), ]

    out = aggregate(shadow_price ~ to, data = sp, FUN = sum)
    names(out) = c("cell", "incoming_dual")

    rc = do.call(rbind, strsplit(as.character(out$cell), "_"))
    out$row = as.integer(rc[, 1])
    out$col = as.integer(rc[, 2])

    out = out[order(-out$incoming_dual, out$row, out$col),
              c("row", "col", "incoming_dual")]
    rownames(out) = NULL
    out
}

settlement_bottlenecks = head(aggregate_bottleneck("shadow_prices_settlement.csv"), 5)
wetland_bottlenecks    = head(aggregate_bottleneck("shadow_prices_wetland.csv"), 5)
```


### Results

We solved the max-flow for both targets:





We solved the max-flow for both targets:







Target



Max flow





Settlement



4.06





Wetland



3.795



Table 1

Settlement has higher max-flow, meaning more fire can reach it - it is harder to protect

Figure 1:  Bottleneck pressure by cell for settlement (left) and wetland (right). Red shading shows the aggregated shadow price for each cell. Bottlenecks concentrate at the target perimeters, not along the current corridor.

### Discussion

The current east-west corridor along row 10 does not coincide with the bottleneck corridors. The shadow prices show fire pressure is highest at rows 6-7, right next to the targets. The corridor is too far south and leaves two cells open on each side, so fire can bypass it.



Settlement is harder to protect than wetland, with a higher max-flow of 4.06 compared to 3.795. This is because the east side of the grid has denser forest vegetation, which gives higher spread probabilities toward the settlement.



This model treats spread probabilities as fixed capacities, so it ignores the stochastic nature of fire spread. It also assumes no wind and a fixed southern ignition front. In practice, wind can shift the direction of fire flow significantly, and ignition can occur anywhere. These limitations are addressed in Steps 2 and 3.



Based on the shadow prices, we would propose clearing cells along column 15 (rows 3-7) to protect the settlement, and along column 7 (rows 3-7) for the wetland. These form barriers at the points where fire pressure is highest. However, this does not account for budget constraints, which we introduce in Step 2.



 

## Step 2

### Model

With a limited budget of B cells, we formulate a MILP to choose which cells to clear as firebreaks to minimise total weighted damage. Settlement cells have weight 10 and wetland cells weight 5. Ignition cells (row 21) and target cells cannot be cleared.



We try two propagation models that share the same budget constraint and objective:



$$\min \sum_{t \in \text{TARGETS}} w_t \cdot d_t \quad \text{s.t.} \quad \sum_{v} \text{clear}_v \leq B
$$


Model A (Binary propagation): An edge(u,v) is active if $p_{\text{burn}}(u,v) > \tau$. If a burning neighbour has an active edge to cell $v$ and $v$ is not cleared, then$v$ must burn:



$$\text{burn}_v \geq \text{burn}_u - \text{clear}_v \quad \forall (u,v) : p_{\text{burn}} > \tau
$$


Every active edge carries fire with certainty, so $\tau$ controls conservatism.

Model B (Continuous reachability): Instead of binary burn, we track a continuous reachability $x_v \in [0,1]$for each cell. Fire propagates multiplicatively with compressed edge weights $$w_{uv} = p_{\text{burn}}(u,v)^{\alpha}$:



$$x_v \geq w_{uv} \cdot x_u - \text{clear}_v \quad \forall (u,v) \in E
$$


The compression$\alpha < 1$ prevents reachability from decaying below solver tolerance over many hops. The objective directly uses reachability at targets:$\min \sum_t w_t \cdot x_t$.



We solve both models at budget levels $B = 5, 8, 10, 13, 16, 20$, sweeping $\tau = 0, 0.05, \ldots, 0.20$ for Model A and $\alpha = 0.05, 0.1, 0.15, 0.2, 0.3$ for Model B.



### Code

Model A — binary propagation MILP (.mod):

```
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

Model B shares the same sets, budget, ignition, and target constraints as Model A.

The differences are: burn is replaced by a continuous reachability x ∈ [0,1], edge weights are compressed as w = p_burn^α to prevent decay over many hops, and propagation is multiplicative rather than binary.

```
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

Both models share a single .run skeleton:

```
option solver gurobi;
option gurobi_options 'mipgap=0.01 timelim=120';
solve;

printf "SUMMARY,%f,%d\n", Total_Damage, solve_result_num;
for {v in CELLS: clear[v] > 0.5} printf "CLEARED,%s\n", v;
for {t in TARGETS: burn[t]  > 0.5} printf "DAMAGED,%s,%d,%f\n", t, weight[t], burn[t];
```

The Model B .run  replaces the burn[t] > 0.5 filter with x[t] > 0.01.



The shared .dat is written from the Step 1 edge frame, targets.csv, and the landscape matrix:

```
#' Write shared AMPL .dat for Step 2 (CELLS, IGNITION, TARGETS, weight, EDGES + p_burn).
#'
#' @param edges_df  from spread_probabilities(): from_row, from_col, to_row, to_col, p_burn
#' @param targets   data frame: row, col, weight
#' @param landscape 21x21 integer matrix (0 = bare, non-zero = vegetated)
#' @param outfile   path for the .dat file
write_step2_dat = function(edges_df, targets, landscape, outfile = "step2.dat") {
    node = function(r, c) paste0(r, "_", c)
    nr   = nrow(landscape)

    veg      = which(landscape != 0, arr.ind = TRUE)
    CELLS    = node(veg[, "row"], veg[, "col"])
    IGNITION = CELLS[veg[, "row"] == nr]
    TGTS     = node(targets$row, targets$col)

    E = data.frame(
        from   = node(edges_df$from_row, edges_df$from_col),
        to     = node(edges_df$to_row,   edges_df$to_col),
        p_burn = edges_df$p_burn
    )
    E = E[E$from %in% CELLS & E$to %in% CELLS, ]

    con = file(outfile, open = "w")
    on.exit(close(con))

    writeLines(c("set CELLS :=",    paste(" ", CELLS,    collapse = ""), ";"), con)
    writeLines(c("set IGNITION :=", paste(" ", IGNITION, collapse = ""), ";"), con)
    writeLines(c("set TARGETS :=",  paste(" ", TGTS,     collapse = ""), ";"), con)

    writeLines("param weight :=", con)
    writeLines(sprintf("  %s  %d", TGTS, targets$weight), con)
    writeLines(";", con)

    writeLines("param: EDGES: p_burn :=", con)
    writeLines(sprintf("  %s  %s  %.8f", E$from, E$to, E$p_burn), con)
    writeLines(";", con)
}

write_step2_dat(edges, targets, landscape, outfile = "step2.dat")
```


The AMPL solve is driven from R via a per-call driver .run so B and the swept parameter are set externally:

```
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
               damage = as.numeric(parts[2]),
               n_cleared = length(cleared),
               cleared_cells = paste(cleared, collapse = ";"),
               status = as.integer(parts[3]))
}
```


Driver loops that produce the sweep results and the refined tau check:

```
A_sweep = NULL
for (b in c(5, 8, 10, 13, 16, 20))
    for (t in c(0, 0.05, 0.10, 0.15, 0.20))
        A_sweep = rbind(A_sweep, run_one("A", b, "tau", t))

A_refined = run_one("A", 13, "tau", 0.105)

B_sweep = NULL
for (b in c(5, 8, 10, 13, 16, 20))
    for (a in c(0.05, 0.10, 0.15, 0.20, 0.30))
        B_sweep = rbind(B_sweep, run_one("B", b, "alpha", a))
```


Marginal value is $(D(B_j) - D(B_{j+1})) / (B_{j+1} - B_j)$within each sweep. Each optimal plan - plus the current corridor (row 10, cols 3-19), the Step 1 bottleneck plan (col 15 rows 3-7 + col 7 rows 3-7), and the no-firebreak baseline - is validated with $K = 1000$ no-wind simulations:

```
#

#' Mean weighted damage + 95% CI for one firebreak plan.
#'
#' @param landscape 21x21 matrix with firebreaks applied (cleared cells -> 0)
#' @param ignition  two-column (row, col) ignition matrix
#' @param targets   data frame: row, col, weight
#' @param K         number of replications
validate_plan = function(landscape, ignition, targets, K = 1000) {
    set.seed(1)
    damages = replicate(K, {
        sim = simulate_fire(landscape, ignition, wind_speed = 0, wind_dir = 0)
        compute_damage(sim$burned, targets)
    })
    mu = mean(damages)
    se = sd(damages) / sqrt(K)
    c(mean = mu, ci_lo = mu - 1.96 * se, ci_hi = mu + 1.96 * se)
}
```




### Results

We solved both models across all budget levels and validated each optimal plan with K  = 1000 fire simulations ( no wind, southern-front ignition, set.seed(1)).



Model A exhibits a sharp threshold effect. In the main sweep, at $\tau \in {0, 0.05, 0.10}$the threshold graph still contains ignition-to-target paths and damage remains 240; at $\tau \in {0.15, 0.20}$ those paths disappear and damage falls to 0. A refined check at $\tau = 0.105$ $ already returns damage 0, placing the transition in $(0.10, 0.105]$. The only useful case is therefore $\tau = 0$ at $B = 13$, Model A places all 13 firebreaks along column 15 and row 7 to protect settlement, sacrificing wetland entirely (MILP damage = 80, simulated mean = 4.8).



Model B produces a smooth budget-damage curve. At$\alpha = 0.1$:













Budget B



MILP damage



Sim. mean(K = 1000)



95% CI





0 (baseline)



—



16.4



[15.1, 17.7]





5



38.2



10.1



[9.1, 11.1]





13



14.1



4.8



[4.3, 5.4] 





16



12.0



3.0



[2.5, 3.4]





20



0



0



[0, 0]





Corridor (17)



—



1.6



[1.1, 2.1]





Per-cell gains are front-loaded:















Budget jump



$$\Delta D$ / cell (sim)$





$0 \to 5$



$	1.26$





$5 \to 13$



$	0.66$





$13 \to 16$



$	0.60$





The Step 1 bottleneck plan (10 cells) improves on baseline but is weaker than Model B at $B = 13$ because diagonal detours remain open; Model A at $B = 13$ , $\tau = 0$ recovers the same settlement-only and settlement-focused geometry.

Figure 2: Left: Model B damage vs budget for differentα velues. Centre: MILP predictions vs simulator mean with 95% CI. Right: firebreak location maps for the current corridor, Model A $(B=13,τ=0),$ and Model B $(B=13,α=0.1)$.





### Discussion

Model A  is effectively all-or-nothing. At $\tau = 0$ every edge is active, so only a complete cut drives $D$ below 240; partial barriers still leave active routes around the firebreak. In the tested budget sweep, the first budget at which such a settlement-covering cut appears is $$ B=13 $$ ,  where Model A places the col-15 + row-7 wall and protects settlement while leaving wetland exposed.



Both deterministic models overestimate damage relative to the simulator ( for example, Model A gives 80 versus simulated mean $4.8$at $B = 13$ ） The MILP treats every active edge as if fire propagates with certainty, so it captures the worst-case superposition of all feasible paths, whereas the simulator realises only one stochastic spread pattern per run.  The deterministic models are therefore used as placement heuristics rather than literal damage predictors: they identify high-pressure corridors, and the resulting geometries are then checked against simulation.


An LP-style shadow price on the budget constraint is not well-defined here because clear is binary, so the integer optimum need not vary smoothly with $B$. For that reason, Table 4 reports adjacent-budget marginal gains $\Delta D / \Delta B$ as an empirical budget-sensitivity proxy instead.



Under Model B, these gains are front-loaded: the simulated reduction per additional cell is much larger over the first few cells than after the main settlement barrier is in place. Model A gives the same qualitative message more sharply, with a threshold jump at $B = 13$. Taken together, the two models suggest that $B = 13$ is the lowest tested budget at which a serious permanent settlement-focused barrier becomes available, while larger budgets mainly extend protection toward wetland.

 

For Model B, the compression exponent $\alpha$is chosen by experimentation: $\alpha = 0.5$ lets reachability decay below solver tolerance over long paths, while $\alpha = 0.1$ preserves useful signal across the full grid. At $B = 13$, the same L-shape is solver-certified optimal for $\alpha \in {0.05, 0.10, 0.15}$, and it remains the incumbent at $\alpha = 0.20$. A different incumbent appears at $\alpha = 0.30$, but that solve terminates with status 402 and is not certified optimal.



We therefore carry the L-shape into Step 3 as the best-supported $B = 13$ permanent candidate on certification and structural-stability grounds, while treating the $\alpha = 0.30$ alternative as unresolved sensitivity rather than as evidence for a different final policy.



An equal-weights sensitivity check at $5{:}5$ keeps the same L-shape as the Gurobi incumbent at both 120 s and 600 s, with the same objective value $14.07$, but both solves still end with status 402 and substantial relative MIP gaps ($25.2\%$ and $14.0\%$). This does not prove that weights are irrelevant; it shows only that, within the tested solver budget, no better balanced geometry displaced the settlement-focused incumbent. The result therefore suggests that the Step 2 $B = 13$scaffold is influenced not only by the 10:5 weighting but also by the underlying south-front geometry.



All Step 2 plans are still calibrated to no wind and a row-21 ignition front, so their robustness under stochastic wind and ignition is deferred to Step 3.



## Step3

### Model

We evaluate the best-supported $B = 13$ permanent plan from Steps 1–2 under stochastic wind and random ignition. This is the L-shape from Step 2 Model B with $\alpha = 0.1$: column 15, rows 1–7, and row 7, columns 16–21. It is solver-certified optimal for $\alpha \in {0.05, 0.10, 0.15}$and remains the incumbent at $\alpha = 0.20$, so we take it into Step 3 as the stress-test candidate rather than as a final policy verdict.



Wind speed is sampled from $V \sim \text{Weibull}(k = 2, A = 7)$ by inverse transform. Wind direction follows a mixture follows a mixture of two von Mises:



$$f_\Theta(\theta) = 0.7 \cdot \frac{e^{3\cos(\theta - 1.5\pi)}}{2\pi\, I_0(3)} + 0.3 \cdot \frac{e^{2\cos(\theta - 0.75\pi)}}{2\pi\, I_0(2)}
$$


The dominant mode ($\mu_1 = 1.5\pi$, westerly) pushes fire toward the settlement; the secondary mode ($\mu_2 = 0.75\pi$, south-easterly) toward the wetland. The von Mises CDF has no closed form, so we use acceptance-rejection with Uniform $(0, 2\pi)$ proposal, envelope $M = \max_\theta f_\Theta(\theta)$, and acceptance condition $U<f 
Θ

 (θ)/M$. The uniform proposal is simple and exact on a bounded support, though it is not efficiency-optimal; a single von Mises with intermediate location and concentration would likely be more efficient.



Each of $N$ scenarios samples $(V, \Theta, \text{ignition})$, runs the simulator, and records damage. The three plans are evaluated on the same sampled batch of wind and ignition scenarios, which controls for exogenous scenario variation, but the simulator's internal random spread draws are not pathwise aligned across plans. We report$\mathbb{E}[D]$ with $95\%$ CI and ${CVaR}_{0.9} = \mathbb{E}[D \mid D \geq q_{0.9}]$ with bootstrap CI. We choose $N=5000$ because it gives adequate precision for the broad comparison at the observed variability; a smaller pilot run gave a similar estimate for the proposed plan.



### Code

Three files implement the Monte Carlo evaluation: step3_samplers.R, step3_monte_carlo.R, and step3_run.R. The functions simulate_fire(), set_firebreaks(), and compute_damage() are from the provided course resource fire-simulator.r, used unchanged.



Samplers for wind direction, wind speed, and ignition (step3_samplers.R)

```
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
#'
#' @param n integer, number of samples to draw.
#' @return list with `samples` (numeric length n), `acceptance_rate` (empirical),
#'   `theoretical_rate` (= 1 / (2 * pi * M)), and `M` (envelope constant).
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

#' Samples wind speed (km/h) by inverse transform of Weibull(shape = 2, scale = 7).
#'
#' @param n integer, number of samples to draw.
#' @return numeric vector of length n.
sample_wind_speed = function(n) {
  u = runif(n)
  7 * (-log(1 - u))^(1 / 2)
}

#' Samples ignition cells uniformly over the vegetated cells of `landscape`.
sample_ignition = function(landscape, n) {
  veg_idx = which(landscape != 0, arr.ind = TRUE)
  chosen  = veg_idx[sample(nrow(veg_idx), size = n, replace = TRUE), , drop = FALSE]
  colnames(chosen) = c("row", "col")
  chosen
}
```

Monte Carlo evaluation and risk measures(step3_monte_carlo.R):

```
#' Evaluates all plans on a pre-sampled batch of (wind_dir, wind_speed, ignition)
#' scenarios. Ignition on a cleared cell yields zero damage (no fire starts).
#'
#' @param landscape  21x21 integer matrix (0 = bare, non-zero = vegetated).
#' @param targets    data frame with columns row, col, weight.
#' @param plans      named list; each element is NULL (no firebreaks) or a
#'   two-column matrix of firebreak cells (row, col).
#' @param wind_dir   numeric vector of length N, wind direction in radians.
#' @param wind_speed numeric vector of length N, wind speed in km/h.
#' @param ignition   N x 2 integer matrix with columns row, col.
#' @return named list (one entry per plan) of data frames with columns
#'   scenario, damage, wind_speed, wind_dir, ign_row, ign_col.
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

#' Computes mean damage with a 95% t-interval, and a bootstrap CI for CVaR_alpha.
#'
#' @param damages numeric vector of per-scenario damage values.
#' @param alpha   CVaR tail level (default 0.9).
#' @param n_boot  bootstrap replicates for the CVaR CI (default 2000).
#' @return data frame with rows mean_damage and CVaR_alpha, columns
#'   measure, estimate, ci_lower, ci_upper.
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

Workflow (step3_run.R). Using the landscape matrix and target table defined earlier from landscape.csv and targets.csv:

```
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

## Monte Carlo evaluation on those same scenarios.
results    = run_mc_evaluation(landscape, targets, plans,
                               wind_dir, wind_speed, ignition)
risk_table = do.call(rbind, lapply(names(results), function(nm) {
  cbind(plan = nm, compute_risk_measures(results[[nm]]$damage))
}))
```

### Results

The sampler achieves $M = 0.464$ , acceptance rate 33.8% (theoretical $1/(2\pi M)=34.3\%$ ). KS test $p=0.73$ (Figure 3)

Figure 3: Wind direction samples vs mixture density (left); wind speed samples vs Weibull(2,7) (right)



Using $N = 5{,}000$ scenarios(set.seed(1) ):













Plan



Cells



$\mathbb{E}[D]$



95% CI of $\mathbb{E}[D]$



$\text{CVaR}_{0.9}$



95% CI of $\text{CVaR}_{0.9}$





No firebreaks



0



13.77



[13.19, 14.35]



58.94



[57.48, 66.21]





Corridor



17



4.09



[3.77, 4.42]



23.83



[22.61, 34.28]





Proposed



13



4.34



[4.08, 4.61] 



26.38



[25.31, 27.50]





Figure 4: Damage distributions. Dashed = mean, solid = $CVaR₀.₉$



Figure 5: Proposed plan: damage vs wind speed by compass direction (left); mean damage by ignition cell (right).

A pilot run at $N=1,000$ gave $E[D]=4.69$versus $4.34$ for the proposed plan, supporting the stability of the broad comparison at this Monte Carlo scale.



### Discussion 

Both plans cut$\mathbb{E}[D]$ from $13.8$ to about $4$, but ${CVaR}_{0.9}\approx25$  reveals substantial tail risk.  The corridor is a strong benchmark: its point estimate is lower than that of the proposed 13-cell permanent plan, and the confidence intervals overlap. Step 3 therefore does not show that the permanent $B=13$ plan outperforms the current corridor; it shows that a more targeted 13-cell placement reaches a similar risk scale with fewer permanent cells.



The worst 10% are driven by ignition near targets which can bypass upstream barriers, and by wind regimes that push fire directly toward them. Westerly wind raises settlement-side pressure, while the south-easterly mode leaves the wetland more exposed.



This also indicates how far the no-wind plan is from a wind-aware one. The current L-shape remains sensible under the dominant westerly mode, but a plan chosen with $(V,Θ)$ known in advance would likely shift some permanent protection toward the wetland side under the south-easterly mode.



$E[D] $ and  ${CVaR}_{0.9}$serve different decision purposes. The mean is diluted by the many zero-damage scenarios, whereas ${CVaR}_{0.9}$ focuses directly on severe failures. If the authority is especially concerned with avoiding high-loss outcomes near the targets,  ${CVaR}_{0.9}$is the more informative complement.



Step 2's MILP predicted damage $14.1$ for this plan, compared with stochastic mean $4.34 $here.  The deterministic model remains useful as a placement heuristic, but Step 3 shows that stochastic evaluation can materially change how candidate plans compare.



The choice of $N$depends on the required precision. Mean-interval width scales roughly as $N^{-1/2}$, whereas tail measures such as ${CVaR}_{0.9}$ require more scenarios because they use only the worst $10\%$ of outcomes.  Here, $N=5000$ is sufficient for a stable broad comparison, but not for claiming a resolved ordering when intervals overlap.





Whether this stress-test verdict supports adopting the 13-cell permanent plan as final policy — in place of the corridor or as part of a two-stage split — is deferred to the Proposal and top-level Discussion, which bring in Step 4 and the authority's operational preferences.



## Step 4

### Model

Step 4 splits the combined budget $B + R = 17$ into two stages. The $B = 13$ permanent cells from Step 3 are fixed before the fire season. After wind $(V, \Theta)$ is revealed but before ignition, up to$R = 4$ additional cells may be cleared reactively. The reactive MILP observes wind only, not the realised ignition $I$; it therefore requires an ignition surrogate（below.）



The reactive cells are chosen by a per-scenario MILP that extends Model B. Permanent cells enter as fixed data (set PERM). To avoid a symbol clash with edge weights, target weights are written $c_t$ (wetland 5, settlement 10). The reachability variable $x_v$ is continuous in$ [0, 1] $— the compressed-propagation formulation from Step 2 Model B, not a binary burn indicator. The only new decisions are binary$z_v$$ for $$  $v\in\mathrm{CELLS}\setminus\mathrm{PERM}$ , with



$$w_{uv} = p_{\mathrm{burn}}(u \to v;\, V, \Theta)^{\alpha}$$


the wind-conditioned edge weight ($α = 0.1$ from Step 2):



$$\min \; \sum_{t \in \mathrm{TARGETS}} c_t\, x_t$$




subject to Model B's propagation / blocking on non-permanent cells, plus the reactive layer:





proxy ignition: $x_u = 1$ for $u \in \mathrm{IGNITION}$



permanent cleared: $x_v = 0$ for $v \in \mathrm{PERM}$



reactive budget: $\sum_{v \in \mathrm{CELLS} \setminus \mathrm{PERM}} z_v \le R$



no reactive on protected sets: $z_v = 0$ for $v \in \mathrm{IGNITION} \cup \mathrm{TARGETS} \cup \mathrm{PERM}$ 



reactive blocking: $x_v \le 1 - z_v$



reactive propagation: $x_v \ge w_{uv}\, x_u - z_v$ for $(u,v)\in\mathrm{EDGES},\;u,v\in/\mathrm{PERM}$



Model B's propagation and blocking on non-permanent cells carry over unchanged and are not repeated. I is unobserved, so the MILP needs a surrogate for$ x_u$. We inherit Step 2's southern-front pattern ($x_u$ = 1 for row 21); this biases the optimisation toward southern geometries (see Q6).



We compare four plans under the same $N = 1{,}000$scenarios sampled with set.seed(1) : (1) no firebreaks, (2)$B = 13$ permanent only, (3) $B + R = 17$ all-permanent,and (4) $B = 13$ permanent $+, R = 4$reactive. Within each scenario, set.seed(1000 + i)  aligns the simulator's starting RNG state across plans, reducing noise in per-scenario comparisons without creating exact pathwise coupling ( different firebreak layouts consume runif at different rates).



### Code

Step 4 reuses Step 3's samplers and compute_risk_measures, and adds three roxygen-documented functions: edges_wind, solve_reactive, and evaluate_all_plans. Only their cores are shown.





Reactive MILP — step4_reactive.mod. Sets (CELLS, IGNITION, TARGETS, EDGES) and param weight {TARGETS} carry over from Step 2 Model B; the additions for Step 4 are shown below.

```
set PERM within CELLS;                   # permanently cleared cells (fixed)
param w {EDGES};                         # p_burn^alpha, pre-computed in R
param R >= 0;                            # reactive clearing budget

var react {CELLS} binary;                # 1 if cell is reactively cleared
var x     {CELLS} >= 0, <= 1;            # continuous fire reachability

check {u in IGNITION}: u not in PERM;
check {t in TARGETS}:  t not in PERM;

minimize Total_Damage: sum {t in TARGETS} weight[t] * x[t];

subject to Reactive_Budget:                         sum {v in CELLS} react[v] <= R;
subject to Ignition_Reach   {u in IGNITION}:        x[u]     = 1;
subject to Ignition_NoReact {u in IGNITION}:        react[u] = 0;
subject to Target_NoReact   {t in TARGETS}:         react[t] = 0;
subject to Perm_Reach       {v in PERM}:            x[v]     = 0;
subject to Perm_NoReact     {v in PERM}:            react[v] = 0;
subject to React_Blocks     {v in CELLS diff PERM}: x[v] <= 1 - react[v];
subject to Propagation {(u,v) in EDGES: u not in PERM and v not in PERM}:
    x[v] >= w[u,v] * x[u] - react[v];
```



edges_wind — wraps spread_probabilities with α-compression:

```
edges_wind = function(landscape, wind_speed, wind_dir, alpha) {
  edges   = spread_probabilities(landscape, wind_speed = wind_speed,
                                 wind_dir = wind_dir)
  edges$w = edges$p_burn ^ alpha
  edges
}
```




solve_reactive — writes .dat + .run, calls AMPL+Gurobi, parses tagged printf output and relmipgap= from Gurobi's log:

```
solve_reactive = function(landscape, targets, perm_cells,
                          wind_speed, wind_dir, alpha, R = 4L,
                          mod_file = "step4_reactive.mod",
                          timelim = 120, mipgap = 0.01, threads = 1L,
                          ampl_bin   = Sys.getenv("AMPL_BIN"),
                          gurobi_bin = Sys.getenv("GUROBI_BIN")) {

  dat = tempfile(fileext = ".dat"); on.exit(unlink(dat), add = TRUE)
  run = tempfile(fileext = ".run"); on.exit(unlink(run), add = TRUE)
  build_step4_dat(landscape = landscape, targets = targets,
                  perm_cells = perm_cells,
                  wind_speed = wind_speed, wind_dir = wind_dir,
                  alpha = alpha, R = R, outfile = dat)

  writeLines(c(
    "reset;",
    sprintf("model %s;", mod_file),
    sprintf("data %s;",  dat),
    sprintf("option solver '%s';", gurobi_bin),
    sprintf("option gurobi_options 'timelim=%d mipgap=%g threads=%d';",
            timelim, mipgap, threads),
    "solve;",
    'printf "SUMMARY,%f,%d\\n", Total_Damage, solve_result_num;',
    'for {v in CELLS: react[v] > 0.5} { printf "REACTED,%s\\n", v; }',
    'for {t in TARGETS: x[t] > 0.01}  { printf "DAMAGED,%s,%d,%f\\n", t, weight[t], x[t]; }'
  ), run)

  out = system2(ampl_bin, run, stdout = TRUE, stderr = TRUE)

  # SUMMARY -> damage, status
  summary_line = grep("^SUMMARY,", out, value = TRUE)
  if (length(summary_line) == 0) {
    return(list(damage = NA_real_, status = NA_integer_, mip_gap = NA_real_,
                react_cells = character(0), n_react = 0L, damaged = NULL))
  }
  parts  = strsplit(summary_line[1], ",")[[1]]
  damage = as.numeric(parts[2]); status = as.integer(parts[3])

  # REACTED -> cell names
  react_cells = sub("^REACTED,", "", grep("^REACTED,", out, value = TRUE))

  # relmipgap= from Gurobi log (last occurrence)
  gap_line = grep("relmipgap=", out, value = TRUE)
  mip_gap  = if (length(gap_line) > 0) {
    as.numeric(sub("relmipgap=", "",
      regmatches(gap_line[length(gap_line)],
                 regexpr("relmipgap=[0-9.eE+-]+", gap_line[length(gap_line)]))))
  } else if (status == 0) 0 else NA_real_

  list(damage = damage, status = status, mip_gap = mip_gap,
       react_cells = react_cells, n_react = length(react_cells))
}
```

evaluate_all_plans — shared-scenario paired evaluator. Output is a long data frame with rows indexed by (scenario, plan). Ignition on a firebreak is detected by land[ir, ic] == 0, short-circuiting the simulator (and, for plan 4, the MILP). Each plan resets the RNG to sim_seed so the four streams start aligned; they drift as firebreak layouts change runif consumption inside simulate_fire.

```
# helper: parse "r_c" node strings to a 2-column (row, col) integer matrix
parse_nodes = function(nodes) {
  if (length(nodes) == 0) return(matrix(integer(0), ncol = 2,
                                        dimnames = list(NULL, c("row","col"))))
  rc = do.call(rbind, lapply(strsplit(nodes, "_"), as.integer))
  colnames(rc) = c("row", "col"); rc
}

evaluate_all_plans = function(landscape, targets, perm_13, perm_17,
                              scenarios, alpha, R = 4L, ...) {

  N          = nrow(scenarios)
  plan_names = c("none", "perm_13", "perm_17", "reactive")
  total_rows = N * length(plan_names)

  land_none = landscape
  land_13   = set_firebreaks(landscape, as.matrix(perm_13[, c("row","col")]))
  land_17   = set_firebreaks(landscape, as.matrix(perm_17[, c("row","col")]))

  out_damage = numeric(total_rows);  out_gap   = rep(NA_real_,    total_rows)
  out_react  = vector("list", total_rows)
  out_plan   = character(total_rows); out_scen = integer(total_rows)

  for (i in seq_len(N)) {
    ir = scenarios$ign_row[i];  ic = scenarios$ign_col[i]
    V  = scenarios$wind_speed[i]; dir = scenarios$wind_dir[i]
    ign      = matrix(c(ir, ic), 1, dimnames = list(NULL, c("row","col")))
    sim_seed = 1000L + i

    # plans 1-3: fixed landscapes — identical structure, pasted once
    for (k in seq_along(list(land_none, land_13, land_17))) {
      land = list(land_none, land_13, land_17)[[k]]
      idx  = (i - 1) * 4L + k
      out_scen[idx] = i; out_plan[idx] = plan_names[k]
      out_damage[idx] = if (land[ir, ic] == 0) 0 else {
        set.seed(sim_seed)
        compute_damage(simulate_fire(land, ign, V, dir)$burned, targets)
      }
    }

    # plan 4: reactive — skip MILP if ignition lands on a permanent cell
    idx = (i - 1) * 4L + 4L
    out_scen[idx] = i; out_plan[idx] = "reactive"
    if (land_13[ir, ic] == 0) {
      out_damage[idx] = 0
    } else {
      sol = tryCatch(
        solve_reactive(landscape, targets,
                       as.matrix(perm_13[, c("row","col")]),
                       V, dir, alpha, R = R, ...),
        error = function(e) NULL)
      if (!is.null(sol) && !is.na(sol$damage)) {
        out_gap[idx]     = sol$mip_gap
        out_react[[idx]] = sol$react_cells
        land4 = set_firebreaks(land_13, parse_nodes(sol$react_cells))
        out_damage[idx] = if (land4[ir, ic] == 0) 0 else {
          set.seed(sim_seed)
          compute_damage(simulate_fire(land4, ign, V, dir)$burned, targets)
        }
      } else {
        out_damage[idx] = NA_real_
      }
    }
  }

  data.frame(scenario_id = out_scen, plan = out_plan,
             damage = out_damage, mip_gap = out_gap,
             react_cells = I(out_react), stringsAsFactors = FALSE)
}
```

Sensitivity runs (R ∈ {0, 2, 5, 6, 8} and (B = 12, R = 5) at N = 200) call evaluate_all_plans with different (perm, R). CVaR0.9CVaR0.9 CIs use compute_risk_measures (2,000 bootstrap replicates).

### Results



We evaluate the four plans on $N=1000$ shared scenarios. $N=1000$ balances precision against runtime: each scenario requires a reactive MILP solve, so Step 4's per-scenario cost is roughly two orders of magnitude larger than Step 3's. At this $N$ the ordering is stable — the no-firebreak and $B = 13$ plans are clearly worse than the two 17-cell plans — but the two 17-cell plans overlap on both mean and $CVaR$ intervals, and $B = 13$ still overlaps them on $CVaR$, so we do not claim a fully resolved ordering.













Plan



$$ \mathbb{E}[D] $(95\%t-CI)$



$$ \text{CVaR}_{0.9} $(95\%bootstrapCI)$





No firebreaks



$14.27 \pm 1.32$



$$64.47$ [$55.10$, $69.25$]$





$B = 13$ permanent



$4.84 \pm 0.60$



$$28.45$ [$24.57$, $30.55$]$





$B + R = 17$ all-permanent



$3.54 \pm 0.54$



$$25.95$ [$18.42$, $28.37$]$





$$B = 13$ + reactive ($R = 4$)$



$3.34 \pm 0.53$



$$19.90$ [$17.75$, $27.10$]$









The scenario-wise gap$D_{\mathrm{perm_13}} - D_{\mathrm{reactive}}$ has mean $+1.49$, reactive strictly reduces damage in $16.8\%$ of scenarios and strictly increases it in $6.9\%$ (plan (2) and (4) shares scenario-level RNG starting states but not exact pathwise realisations.) Against perm_17, reactive is worse in $2.4\%$. Of the $1000$ scenarios, $31$ had ignition on a $B = 13$ permanent cell (zero damage, MILP bypassed); of the remaining $969$ solves, Gurobi reached the$1\%$target on $960$($99.1\%$), with median final gap $0.97\%$.

Figure 6 gives a partial answer to where reactive effort matters most. Ignition location shows the clearest signal: for ignitions in row 1-7 （near the targets) the mean gap is $0.37$ with $P(\mathrm{gap}>0)=3.2\%$ ; in rows 8-14 the mean gap rises to $2.68$ ($P(\mathrm{gap}>0)=22.3\%$) ; in rows 15-21 （near the MILP's southern-front proxy) the mean gap is$1.38$ . Wind-speed terciles do not show a monotone trend (mean gaps $1.87 / 0.92 / 1.70$ from low to high), so the evidence does not support the claim that reactive matters most at extreme wind speeds. Against perm_17 specifically, the mean damage difference is essentially zero under westerly wind                               ( $\Delta = 0.03$, $n = 714$) and larger under easterly wind ($\Delta = 0.65$, $n = 286$) . These stratifications primarily reflect where the southern-front proxy aligns or misaligns with the realised fire geometry; they do not yet isolate a pure "marginal value of wind information" signal.

Figure 6: (a) Distribution of the scenario-wise reactive effort gap D_perm_13 − D_reactive with mean (solid) and zero (dashed). (b) Mean gap split by ignition row: reactive improves most when ignition lies in the middle rows (8–14), less at either 



Reactive cell frequency (Figure 7). Four cells dominate: $(11, 4)$, $(11, 5)$, $(11, 6)$ are each selected in $96\%$ of scenarios and $(15,5)$ in $73.7\%$— exactly the four cells the no-wind MILP at $B = 17$ adds on top of $B = 13$ . Of the $3876$ total reactive selections, $99.6\%$ are closer (Chebyshev distance) to the wetland than to the settlement; $0\%$ are closer to the settlement.



Figure 7: Spatial frequency of reactive cell selection across N = 1,000 scenarios. Colour intensity = number of scenarios selecting that cell. Black outlines mark the B = 13 permanent plan; red and blue outlines mark the settlement and wetland targets. Reactive effort concentrates on the four cells the no-wind B = 17 MILP adds to B = 13.



R sensitivity (Figure 8).  All points use $N = 200$ for consistency. At  $R=0$ plan (4) reduces to  perm_13 ($E=3.70$ \text{, } $CVaR=17.24$); $R=2$ lowers the mean to $3.47$ but $CVaR_{0.9}\left(24.05\right)$is anomalous .$R = 4$ gives $E=3.20$, $CVaR=18.17$. ($189 / 193$ optimal, median gap $0.97\%$ )$R = 5$: $2.57/17.69(131/193,1.00\%)$ ; $R = 6$ : $1.25/10.87(5/193,27.3\%)$; $R = 8$ :$1.47/10.93(143/193,0.98\%)$ $$R = 6$ is a solver outlier — the combinatorial choice at that budget is especially ill-conditioned for Gurobi under the $120$ s limit — but $R = 5$ and $R = 8$ both finish near$1\%$, so the downward trend from $R = 4$ onward is not an artefact of one bad budget. The $R = 6 → 8$ non-monotonicity is within $N = 200$ sampling noise. Overall, diminishing returns around $R = 6–8$; the exact elbow is not pinned down.

Figure 8: R sensitivity of plan (4), all at N = 200. (a) Mean damage and (b) $\text{CVaR}_{0.9}$ as functions of the reactive budget $R \in {0, 2, 4, 5, 6, 8}$, with perm_13 and perm_17 as horizontal references. The R = 6 point is circled because its MIP gap is large and the result is not certified optimal.





### Discussion

Q1: A fully optimal approach is a two-stage stochastic program with non- anticipative recourse, where the reactive policy$z(V, \Theta)$ depends on the revealed wind but not on the unobserved ignition$I$ : 





$$\min_{\mathrm{perm},\ z(\cdot)}\ \mathbb{E}_{V,\Theta}\bigl[\ \mathbb{E}_{I}\bigl[\ D(\mathrm{perm},\ z(V,\Theta);\ V, \Theta, I)\ \bigr]\ \bigr].$$


This requires scenario-based decomposition over the joint $(V, \Theta, I)$ distribution,  Since Step 3 samples ignition uniformly over the 439 vegetated cells and $(V,\Theta)$ is continuous, the exact formulation is computationally at a feasible budget. Our Step 4 surrogate replaces the inner expectation with a southern-front proxy and solves one MILP per sampled $(V, \Theta)$.



Q2: Wind information weakly enlarges the reactive policy class, so one would expect plan (4) to match or beat plan (3) under an exact two-stage model; the surrogate does not guarantee exact sample-level dominance. The observed means are consistent ($3.34$ vs $3.54$ ; $CVaR₀.₉ $19.90$$ vs $25.95$), driven by the easterly-wind regime $($\Delta = 0.65$, $n = 286$)$; under westerly wind the  plans are indistinguishable $($\Delta = 0.03$, $n = 714$)$ . The $2.4\%$ of scenarios where reactive exceeds perm_17 are consistent with the ignition surrogate (Q6) and residual simulator -RNG drift rather than a systematic modelling failure.


Q3: The$99.1\%$ rate above implies the optimality gap is unlikely to drive reactive cell choices in our sample. Separately, the $B = 17$ all-permanent MILP reported a $25\%$ gap at the 120s limit, but a re-solve at 1800s confirmed the same incumbent is globally optimal — a slow LP lower bound, not a bad incumbent.



Q4: The reactive MILP concentrates its four cells on $(11,4), (11,5), (11,6), (15,5)$ in most scenarios, matching the no-wind MILP's $B = 17$ additions to $B = 13$ . Wind induces only small secondary shifts. This suggests the permanent plan leaves one specific structural gap on the wetland side that is almost always worth closing, rather than the reactive stage exploiting fine-grained wind-specific adjustments.



Q5: The settlement weights (10) is twice the wetland weight (5), yet $99.6\%$ of reactive selections are closer to the wetland. This is consistent with the permanent plan already largely saturating settlement protection（ the$B = 13$ L-shape blocks the main northward corridor ), so marginal reduction from additional settlement-side cells is small and reactive effort is redirected toward the wetland. Under this reading, weight asymmetry is absorbed by the permanent plan rather than expressed in the reactive layer.



Q6: The reactive MILP is conditioned on $(V, \Theta)$ only;  $I$ is unobserved, so a surrogate for $x_u$ is required.  We use Step 2's southern-front pattern ($x_u = 1$ for row 21).  The stratified gaps in Results match the expected degradation as ignition moves away from the proxy front: in the northern third reactive rarely helps  $($ P(\text{gap}>0)=3.2 $\%)$ ; middle rows give the largest mean gap ($2.68$). A natural alternative — weighting ignition cells by a prior over $I$ — would be more faithful to a stochastic-program formulation; it changes only the MILP's input, leaving outer MC and bootstrap $CVaR$ CIs unchanged. We keep the current surrogate to match Step 2 and keep the MILP small, but this is a modelling choice, not a principled optimum.



Q7: Re-solving Step2 Model B at $B = 12$ gives  perm_12 = perm_13  minus cell$(1, 15)$.  $(B=12,R=5)$at $N = 200$ gives $\mathbb{E} = 3.15$$[2.07,4.23]$and ${CVaR}_{0.9} = 18.45$$[13.45, 28.75]$ , compared with $3.20$ $[2.12, 4.28]$and $18.17$$[13.89, 28.50]$ for the the matched-$N$ check with $(B = 13, R = 4)$ split,Both mean and $CVaR$ intervals overlap substantially, so trading one permanent cell for one reactive cell does not meaningfully change risk outcomes under our surrogate within this range.

## Discussion

How the mathematics translates into into the proposal

The final proposal is a two-stage plan with a 13-cell permanent L-shape and a reactive budget of 4-cells, matching the corridor's 17-cell count.



Permanent placement. Step 1's max-flow shadow prices peak at rows 6-7 next to the targets, so the row-10 corridor is structurally mis-located. Step 2 converts this into the 13-cell L-shape via Model B $\alpha = 0.1$, which Steps 3-4 carry forward as the permanent layer at $B = 13$. 



Why recourse enters. Step 3 rules out permanent-only at $B = 13$ as a final recommendation: it gives the permanent 13-cell plan $\mathbb{E}[D] = 4.34$ and ${CVaR}_{0.9} = 26.38$ , worse than the current corridor on both metrics despite the smaller footprint. Step 4 adds reactive cells conditioned on the revealed wind; the two-stage plan's point estimates sit below the permanent-only 17-cell plan's , but $95\%$ intervals overlap on both mean and tail, so Step 4 supports preference for the two-stage plan among the 17-cell alternatives, not proven dominance.



Overprediction  vs the simulator: Both deterministic formulations over-estimate absolute damage — Model A reports 80 at $B = 13$ against  simulator mean of 4.8 under no wind — because they represent spread paths deterministically inside the optimisation model, while the simulator realises each propagation attempt stochastically with probability $p_{\text{burn}} < 1$ and many paths die out. the deterministic MILPs are therefore used as placement tools, not as expected- damage predictors.



### Questions for further work

Q1: Is the recommendation actually better than the current corridor?

No direct within-regime head-to-head has been run. Step 3 and Step 4 use different Monte Carlo designs, so their point estimates cannot be differenced across tables. Resolving this needs a single common experiment — the corridor under the Step 4 desig, or the two-stage plan under Step 3's.



Q2. Is the 17-cell total budget justified?

Step 2 shows diminishing returns past $$ B=13 $$ , which justified 13 as the permanent stopping point; the 17-cell cap itself is not separately optimised and matched the current corridor for same-footprint compatison. A full $(B, R)$ sweep beyond the $(13, 4)$ and $(12, 5)$ pairs tested would be needed to claim optimality.



Q3. How sensitive are conclusions to the wind model?

The wind distribution is stylised, and the reactive stage assumes that the wind observed by the MILP matches the simulator wind; the wind-to-$p_{\text{burn}}$ spread mapping is inherited unchanged. Robustness under any of these remains open.



Q4. How much does the southern-front ignition proxy limit Step 4?

The reactive MILP treats all row-21 cells as fully reachable because ignition is unobserved decision time. Step 4 Figure 6 shows the reactive-vs-permanent-13 gap is smallest for northern-third ignitions, consistent with the proxy being least accurate there.



Q5. What would resolve permanent-17 vs two-stage? 

At $$ N=1000 $$ the two 17-cell plan's intervals overlap on both metrics, A paired comparison on scenario-wise differences (using Step 4's shared- scenario design) or a lager $N$ would be needed before upgrading preference to dominance. The $R$-sweep past $R = 4$ is suggestive, not decision-grade, because those runs use $N = 200$ with MIP quality degrading for $R \geq 6$.



### Summary

Current evidence supports the two-stage plan as the best-supported evaluated alternative; it does not prove superiority over either the current corridor or the permanent-only 17-cell plan.





## AI Declaration

AI use. I use Claude to help draft initial AMPL templates for the Step 1 LP and the Step2/Step 4 MILPs, which i then rewrote, checked ,and aligned with my data interface and constraint logic, and to review my R code for indexing, RNG handling across scenarios, and edge cases, plus minor plot and LaTex formatting.



What was mine. Modelling choices (network max-flow; Model A binary vs Model B continuous-reachability MILP; two-stage stochastic recourse); the compression exponent $\alpha = 0.1$, chosen after$\alpha = 0.5$ underflowed across long propagation paths; the beyond-course sampling design (mixture of two von Mises component with Bessel $I_0$ normalisation, envelope $M$ for acceptance-rejection, acceptance-rate and KS diagnostics); and all result interpretation. Reactive cap $R = 4$ follows the Step 4 brief. Permanent budget $B = 13$ is the Step 3-4 task setting; Step 2 separately evaluates a range of $B$ values and my analysis supports carrying $B = 13$ forward at the marginal-gain elbow. All prose in the submission was written be me; AI edits were limited to delimiters and grammar.



Reflection. The discipline I keep from this workflow is explicit set scope and constraint activation conditions when reviewing AMPL. It did not change any modelling decision.
