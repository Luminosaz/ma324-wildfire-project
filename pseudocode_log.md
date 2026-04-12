# pseudocode / formulation notes

my working notes before asking Claude for help implementing stuff

---

## step 1 — max-flow

model the grid as a network, fire flows from south (row 21) toward targets

nodes: all vegetated cells + source S + sink T
edges: adjacent cells (8 directions), capacity = p_burn from simulator
S connects to all row 21 cells (cap = big number)
target cells connect to T (cap = big number)

standard max-flow with return arc T -> S

```
max flow(T,S)

for all nodes v:
  inflow = outflow  (conservation)

for all edges (i,j):
  flow(i,j) <= cap(i,j)  (capacity — need this as named constraint for shadow prices)
```

shadow prices tell me which edges are bottlenecks
to figure out which CELLS matter: sum up shadow prices of edges going into each cell
run it twice, once for settlement once for wetland

---

## step 2 — MILP firebreak placement

now I have a budget of B cells I can clear. goal: pick which B cells to turn into firebreaks so that total weighted damage to targets is minimised. settlement weight=10, wetland weight=5.

same setup as step 1: no wind, fire starts from entire row 21.

the key challenge is encoding "fire spread" as constraints in an optimisation model. the brief says try two approaches.

both models share:
- clear[v] binary — 1 if I place a firebreak at cell v
- can't clear ignition cells (row 21) or target cells
- budget: total cleared <= B
- if a cell is cleared it can't burn

```
min sum over targets t: weight[t] * burn[t]

sum clear[v] <= B
burn[v] = 1 for ignition cells
burn[v] <= 1 - clear[v] for clearable cells
```

### model A: binary propagation

idea: pick a threshold τ. if p_burn(u→v) > τ, that edge is "active" — meaning fire WILL spread if u is burning and v isn't cleared. edges below τ are ignored.

this makes it worst-case / deterministic: every active edge is treated as certain fire spread.

```
for each active edge (u,v) where p_burn > τ:
  burn[v] >= burn[u] - clear[v]

if u burns and v not cleared → burn[u]=1, clear[v]=0 → burn[v] >= 1
if v is cleared → clear[v]=1 → burn[v] >= burn[u] - 1 → doesn't force anything
```

τ controls how conservative:
- τ=0 means ALL edges active, even very weak ones → super conservative, needs lots of firebreaks
- τ=0.2 ignores anything below 20% → optimistic, less firebreaks needed
- need to sweep τ and see which gives predictions closest to simulator

UPDATE after first sweep:
- τ=0 to 0.1: fire propagates deterministically everywhere, even B=5/8/10 can't block
  anything — model says "don't bother clearing" (n_cleared=0, damage=240)
  only B=13 works because it can cut a full column+row to isolate settlement
- τ=0.15 and 0.2: too many edges removed, fire can't reach targets at all, damage=0
  even with B=0 nothing burns
- problem: no τ value gives a gradual budget-vs-damage curve
- need finer sweep between 0.10 and 0.15 where the transition happens
- revised τ sweep: 0, 0.05, 0.10, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15, 0.20

UPDATE after refined sweep:
- transition is between τ=0.10 and τ=0.105, but it's a cliff not a slope
  τ<=0.10: damage=240 (all targets burn), τ>=0.105: damage=0 (no fire reaches targets)
- this is inherent to binary propagation: it's all-or-nothing
  either fire can reach targets deterministically → damage=max
  or fire path is broken → damage=0, no firebreaks needed
- only τ<=0.1 with B>=13 gives an interesting result (protect settlement, sacrifice wetland)
- conclusion: Model A is a poor tool for budget sensitivity analysis
  but it IS useful for identifying the minimum cut (B=13 at τ=0 → settlement protection)
  this contrast with Model B's smooth curve is a key discussion point

for shadow price: this is a MILP so no LP dual directly. but I can compute marginal value as damage(B) - damage(B+1) by solving at consecutive budgets.

### model B: continuous reachability

instead of binary burn, track a continuous "reachability" x[v] in [0,1] for each cell. x represents how much fire can reach v through the network.

fire propagates multiplicatively: if u has reachability x[u] and edge weight w(u,v), then v gets at least w * x[u] reachability (unless cleared).

problem: raw p_burn values multiply along paths and decay super fast (like 0.25^10 ≈ 10^-6), solver can't distinguish from 0. so compress: w(u,v) = p_burn^α with α < 1.

```
x[v] >= p_burn(u,v)^α * x[u] - clear[v]   for each edge (u,v)
x[v] = 1 for ignition cells
x[v] <= 1 - clear[v] for clearable cells
x[v] >= 0

min sum over targets t: weight[t] * x[t]
sum clear[v] <= B, clear binary
```

nice thing: don't even need a binary burn variable. x[t] IS the damage measure (continuous, between 0 and 1). objective is weighted sum of reachability at targets.

this is MILP because clear is binary, but x is continuous — should solve faster than if everything were binary.

α controls compression:
- α=1 → raw probabilities, might be too small for solver
- α=0.5 → square root, spreads values out
- α=0.3 → heavy compression, everything looks more reachable → conservative

UPDATE after testing model B with gurobi:
- α=0.5 doesn't work — ran B=0 (no firebreaks) and damage was only 0.07
  reachability at targets all < 0.001, model has no discriminating power
  reason: 15+ hops from row 21 to targets, even sqrt(0.45)^15 ≈ 0.003
- tried α=0.1 — B=0 damage=43.5, reachability 0.12-0.24, actually useful
- need α small enough that reachability doesn't vanish before reaching targets
- revised sweep: α = 0.05, 0.1, 0.15, 0.2, 0.3

### what I need to build

```
R scripts:
  step2_build_dat.R — read landscape/targets/edges, write .dat for AMPL
  step2_solve.R — loop over B values, call AMPL for both models, collect results
  step2_simulate.R — run fire simulator to validate model predictions
  step2_plot.R — comparison figures

AMPL files:
  step2_modelA.mod — binary propagation MILP
  step2_modelB.mod — continuous reachability MILP
  step2.run — script to load data, set params, solve, output results

B values to try: 5, 8, 10, 13, 16, 20
τ values: 0, 0.05, 0.10, 0.15, 0.20
α values: 0.05, 0.1, 0.15, 0.2, 0.3
```

---

## step 3 — monte carlo

### AR sampler for wind direction

(todo)

### MC evaluation

(todo)

---

## step 4 — reactive

(todo)
