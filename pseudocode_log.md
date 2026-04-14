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

## step 3 — monte carlo under uncertain wind

step 2 assumed no wind and fire from south edge. step 3 drops that — sample random wind + random ignition point and see how the plan actually holds up.

### distributions (from project brief)

wind speed V ~ Weibull(shape=2, scale=7)
  - easy, R has rweibull() built in
  - mean ~6.2 km/h, occasionally up to 15-20

wind direction θ ~ mixture of two von Mises on [0, 2π):
  - 70% from west (μ₁ = 1.5π, κ₁ = 3) — pushes fire east toward settlement
  - 30% from south-east (μ₂ = 0.75π, κ₂ = 2) — pushes fire north-west toward wetland

  f(θ) = 0.7 * vonmises(θ; μ₁, κ₁) + 0.3 * vonmises(θ; μ₂, κ₂)

  where vonmises(θ; μ, κ) = exp(κ * cos(θ - μ)) / (2π * I₀(κ))
  I₀ = modified bessel function of first kind, order 0

ignition: uniform random vegetated cell

### acceptance-rejection sampler for wind direction

can't use inverse transform because von Mises CDF has no closed form.
acceptance-rejection with uniform proposal on [0, 2π) is the simplest option.

```
target f(θ) = mixture density as above
proposal g(θ) = 1/(2π)   (uniform on [0, 2π))

M = max of f(θ) over [0, 2π)
  find M numerically: evaluate f on dense grid, take max

to sample one θ:
  repeat:
    draw θ ~ Uniform(0, 2π)
    draw u ~ Uniform(0, 1)
    if u < f(θ) / (M * g(θ)):
      accept θ
      break

acceptance rate = 1/M * (1/(2π)) ... wait no.
actually: acceptance prob = integral of f / (M * 2π * (1/2π)) = 1/M
hmm let me think again.

g(θ) = 1/(2π), so M*g(θ) = M/(2π)
we need M*g(θ) >= f(θ) for all θ
f(θ) <= M by definition of M, and g(θ) = 1/(2π)
so M*g(θ) = M/(2π) — but this needs to be >= f(θ) = at most M
that only works if M/(2π) >= M, i.e. 1/(2π) >= 1 — FALSE.

ok I need to use the correct envelope. the envelope is c*g(θ) where c = max(f/g).
since g = 1/(2π), the ratio f/g = 2π*f(θ), so c = 2π*M where M = max f(θ).
acceptance rate = 1/c = 1/(2π*M).

M ≈ max of the mixture density. the von Mises with κ=3 peaks at exp(3)/(2π*I₀(3)) ≈ 20.09/(2π*4.88) ≈ 0.655. so f peaks around 0.7*0.655 ≈ 0.46 or so.
c = 2π * 0.46 ≈ 2.87
acceptance rate ≈ 1/2.87 ≈ 35%

that's decent enough for uniform proposal. could get higher with a von Mises proposal but uniform is simpler to justify.
```

```
sample_wind_dir(n):
  # compute M = max f(θ) on fine grid
  grid = seq(0, 2π, length=10000)
  M = max(f(grid))
  c = 2π * M  # envelope constant for uniform proposal

  samples = empty vector
  n_proposed = 0
  while length(samples) < n:
    θ = runif(1, 0, 2π)
    u = runif(1, 0, 1)
    n_proposed += 1
    if u < f(θ) / (c * (1/(2π))):   # simplifies to u < f(θ)/M
      append θ to samples

  acceptance_rate = n / n_proposed
  return samples
```

wait, let me simplify. the acceptance condition is u < f(θ) / (c * g(θ)).
c = 2π*M, g(θ) = 1/(2π).
c * g(θ) = 2π*M * 1/(2π) = M.
so: accept if u < f(θ) / M.

that's clean.

### MC evaluation loop

```
run_mc_evaluation(landscape, firebreak_cells, targets, N):
  set.seed(1)
  
  # apply firebreaks
  ls_mod = set_firebreaks(landscape, firebreak_cells)
  
  # pre-compute vegetated cells for ignition sampling
  veg_cells = which(landscape != 0, arr.ind = TRUE)
  
  damages = numeric(N)
  wind_speeds = numeric(N)
  wind_dirs = numeric(N)
  ignition_rows = numeric(N)
  ignition_cols = numeric(N)
  
  for i in 1:N:
    V = rweibull(1, shape=2, scale=7)
    θ = sample_wind_dir(1)
    idx = sample(nrow(veg_cells), 1)
    ign = veg_cells[idx, ]
    
    result = simulate_fire(ls_mod, ign, wind_speed=V, wind_dir=θ)
    damages[i] = compute_damage(result$burned, targets)
    wind_speeds[i] = V
    wind_dirs[i] = θ
    ignition_rows[i] = ign[1]
    ignition_cols[i] = ign[2]
  
  return data frame of (damage, V, θ, ign_row, ign_col)
```

### risk measures

```
given damages vector of length N:

E[damage] = mean(damages)
95% CI for E[damage] = mean ± 1.96 * sd(damages) / sqrt(N)

CVaR_0.9 = mean of worst 10% of damages
  sorted = sort(damages)
  tail_start = ceiling(0.9 * N) + 1
  CVaR = mean(sorted[tail_start:N])

95% CI for CVaR: bootstrap
  for b in 1:1000:
    resample = sample(damages, N, replace=TRUE)
    CVaR_b = compute CVaR on resample
  CI = quantile(CVaR_boots, c(0.025, 0.975))
```

### plans to compare

1. no firebreaks (baseline)
2. corridor: row 10, cols 3-19 (17 cells)
3. proposed: Model B, B=13, α=0.1 optimal cells from step 2

### what to look for in discussion

- which wind directions hurt the proposed plan most?
  probably east wind (pushes fire away from where we blocked)
- worst 10%: look at wind dir + ignition location combinations
- corridor (17 cells) vs proposed (13 cells): does broader coverage help in tail risk?
- how far from wind-aware optimal? can't easily compute this, but can argue based on damage distribution shape
- E[damage] vs CVaR: which matters more depends on whether we care about average or worst case
- N justification: show CI width shrinking with N, pick where it stabilises

### sampler validation results (set.seed(1), n=10000)

wind direction AR sampler:
  M = 0.4636 (density peak at θ ≈ 1.5π, western mode)
  acceptance rate = 0.3403, theoretical = 0.3433 (ratio 0.9912)
  KS test p = 0.95 → cannot reject correct distribution
  density integrates to 1 ✓

wind speed Weibull(2,7):
  sample mean = 6.20 (theoretical 6.20)
  sample var = 10.56 (theoretical 10.52)
  KS test p = 0.84 ✓

ignition uniform over vegetated cells:
  439 vegetated cells, 2 bare
  all samples land on vegetated cells ✓
  chi-sq uniformity p = 0.97 ✓

### MC evaluation results (set.seed(1))

N=1000 pilot:
  baseline  E[damage]=13.88 [12.57, 15.18]  CVaR₀.₉=59.73 [56.41, 69.95]
  corridor  E[damage]= 4.10 [ 3.33,  4.86]  CVaR₀.₉=25.79 [22.57, 38.70]
  proposed  E[damage]= 4.69 [ 4.11,  5.27]  CVaR₀.₉=24.46 [22.26, 26.80]

N=5000:
  baseline  E[damage]=13.77 [13.19, 14.35]  CVaR₀.₉=58.94 [57.48, 66.21]
  corridor  E[damage]= 4.10 [ 3.77,  4.42]  CVaR₀.₉=23.83 [22.61, 34.28]
  proposed  E[damage]= 4.35 [ 4.08,  4.61]  CVaR₀.₉=26.38 [25.31, 27.50]

observations:
  - N=1000 and N=5000 estimates consistent → N=1000 already sufficient
  - CI half-width for mean: ~0.8 at N=1000, ~0.3 at N=5000
  - corridor vs proposed E[damage] CIs overlap → difference not significant
  - proposed uses 13 cells vs corridor's 17 → better per-cell efficiency
  - CVaR ≈ 6× mean → heavy tail risk, worst 10% scenarios cause ~60 damage
  - ignition on firebreak cell → damage=0 (handled in code)

---

## step 4 — reactive

two-stage decision. B=13 cells fixed before season (from step 3), plus R=4
reactive cells chosen after wind (V, theta) is revealed but before ignition.

what the reactive stage knows:
  - yes: V, theta (the realised wind for this scenario)
  - no: the actual ignition location (it hasn't happened yet)
  - no: the stochastic outcome of the fire
so the reactive MILP must plan as if ignition could be anywhere consistent
with the prior. brief tells us to use the southern-front proxy — fire from
row 21 — same as step 2. this is conservative: it assumes worst-case
ignition pattern. NOT per-scenario true ignition, that would be cheating.

### the reactive MILP (model B extended)

basically step 2 model B with:
  - permanent cells (13) hard-wired as param, not decision variables
  - new binary var z[v] for reactive cells
  - edge weights w use wind-aware p_burn

```
param perm_clear {CELLS} binary, default 0  # 1 on the 13 permanent cells
var   z          {CELLS} binary              # reactive decision
var   x          {CELLS} in [0,1]            # reachability (continuous)

param p_burn {EDGES}       # recomputed for this scenario's (V, theta)
param alpha  = 0.1         # same as step 2 model B final
param R      = 4
param w {(u,v) in EDGES} = p_burn[u,v] ^ alpha

min sum_{t in TARGETS} weight[t] * x[t]

s.t.
  Reactive_Budget:  sum_v z[v] <= R
  No_Reactive_On_Perm {v}: z[v] <= 1 - perm_clear[v]    # avoid double-clearing
  Ignition_NoClear {u in IGNITION}: perm_clear[u] = 0   # enforced via data
  Ignition_Reach   {u in IGNITION}: x[u] = 1
  Target_NoClear   {t in TARGETS}:
      perm_clear[t] = 0 and z[t] <= 0                   # targets never cleared
  Clear_Blocks {v in CELLS}:
      x[v] <= 1 - perm_clear[v] - z[v]
  Propagation {(u,v) in EDGES}:
      x[v] >= w[u,v] * x[u] - perm_clear[v] - z[v]
```

sanity: if all wind is zero and perm_clear = the step-2 no-wind B=13 plan,
z=empty should be an optimal reactive solution (no gain from R=4 since
perm already blocks). this is a good test before trusting the windy runs.

### wind edges — just wrap the simulator

fire-simulator.r already has `spread_probabilities(landscape, V, theta)`.
returns a data.frame with (from_row, from_col, to_row, to_col, p_burn). so
for step 4, the wind-dependent edges don't need a fresh implementation:

```
edges_wind(landscape, V, theta, alpha):
  edges = spread_probabilities(landscape, V, theta)
  edges$w = edges$p_burn ^ alpha
  return edges
```

golden test before trusting it:
  - V=0 case: edges_wind(L, 0, *, alpha) should match step 2's edges
    (which were built with wind=0). compare row-by-row, tol 1e-10.
  - directionality: for V>0, theta=pi/2 (wind from east) the edges
    pointing west should have higher w than edges pointing east. hand
    check 2-3 edges.

reason to trust simulator's function: step 2 already relied on
`spread_probabilities(L, 0, 0)` to build the no-wind edges, so reusing it
with (V, theta) != 0 is the only way to stay consistent with the sim.

### single-scenario reactive solve — the driver

```
solve_reactive(perm_cells, V, theta, alpha, R, time_limit, mip_gap):
  edges = edges_wind(landscape, V, theta, alpha)
  write perm_clear.dat with 1's on perm_cells, 0 elsewhere
  write edges.dat with w column
  run AMPL: step4_reactive.mod + perm_clear.dat + edges.dat
    option solver gurobi
    option gurobi_options 'timelim=<time_limit> mipgap=<mip_gap>'
    solve
  parse output:
    z_cells = cells where z[v] = 1
    obj = _obj
    solve_time = _solve_elapsed_time
    mip_gap = gurobi's final reported gap
  return list(z_cells, obj, solve_time, mip_gap)
```

### batch evaluation — 4 plans, same scenarios, CRN

```
evaluate_all_plans(scenarios, perm_13, perm_17, alpha, R, solver_opts):
  # scenarios is a data.frame with columns (ign_row, ign_col, V, theta)
  # perm_13 and perm_17 are fixed cell sets

  rows = empty list
  for i in 1..N:
    ign = (scenarios$ign_row[i], scenarios$ign_col[i])
    V   = scenarios$V[i]
    th  = scenarios$theta[i]

    # plan 1: no firebreaks
    sim1 = simulate_fire(landscape, ign, V, th)
    d1   = compute_damage(sim1$burned, targets)

    # plan 2: permanent_13 only
    L2   = set_firebreaks(landscape, perm_13)
    sim2 = simulate_fire(L2, ign, V, th)
    d2   = compute_damage(sim2$burned, targets)

    # plan 3: permanent_17 only
    L3   = set_firebreaks(landscape, perm_17)
    sim3 = simulate_fire(L3, ign, V, th)
    d3   = compute_damage(sim3$burned, targets)

    # plan 4: permanent_13 + reactive z*(V, theta)
    res  = solve_reactive(perm_13, V, th, alpha, R, ...)
    L4   = set_firebreaks(landscape, rbind(perm_13, res$z_cells))
    sim4 = simulate_fire(L4, ign, V, th)
    d4   = compute_damage(sim4$burned, targets)

    rows.append (scenario_id=i, plan=1..4, damage=d1..d4,
                 reactive_cells (NA for 1-3, else res$z_cells),
                 solve_time, mip_gap)
  return tidy data.frame
```

critical: only ONE set of scenarios for all 4 plans — same (ign, V, th)
triplet across plans (CRN). also: same RNG state for simulate_fire within a
scenario across plans (wrap each call in set.seed(324 + i) or similar) —
otherwise the simulator's internal runif noise differs and the comparison
isn't apples-to-apples. will verify in phase D.

### CRN with step 3

step 3 already ran N=5000 scenarios with set.seed(324). for step 4:
  - either persist step 3's (ign, V, theta) to a .rds and load
  - or in step 4's script: set.seed(324); run step 3's sampler with N; use
    those triplets. will give the same triplets IF the sampler code and
    RNG calls are identical. safer to persist, so plan-2 damage in step 4
    must match plan-2 damage in step 3 exactly.
need to check: does step3_run.R dump the scenarios? if not, dump them as
step 3 finalisation before starting step 4 so numbers line up.

### plan (3): pre-season at B+R = 17

same step 2 model B (alpha=0.1, no wind), solved at B=17 instead of 13.
gives a 17-cell permanent-only plan. then evaluate with the same 4-plan
machinery above. implementation: just call step 2's solve script with B=17,
read the 17 cells.

### R sensitivity and (B=12, R=5)

R sweep: R in {0, 2, 4, 6, 8}, same N=500 pilot (R=4 bumped to N=1000).
for R=0 plan 4 degenerates to plan 2 — automatic check.

B=12, R=5: re-solve step 2 model B at B=12 (new 12-cell permanent), then
run the whole 4-plan evaluation with R=5. cost: ~1 extra batch.

### risk measures

reuse step 3's compute_risk_measures: E[damage], CVaR_0.9, t-CI for mean,
bootstrap CI for CVaR. do not rewrite.

### phase order (from the plan file, summarised)

B: get edges_wind + golden test working. nothing downstream works without this.
C: single-scenario solve_reactive working + hand-verify z is sensible.
D: full 4-plan batch with CRN.
E: R sweep and B=12/R=5.
F: figures + report.

open questions (still need to resolve before phase B):
  - does step 3 persist the 13 cells? if not, dump from step2_results.
  - does step 3 persist the N=5000 scenarios? if not, dump before touching
    step 4 so CRN is preserved.
  - plan 3: B=17 model B at alpha=0.1 — run step 2's solver, same as
    before just with B=17. needs doing in phase D.
