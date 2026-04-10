# MA324 Wildfire Firebreak Planning

This is my code repo for the MA324 summative project (LSE, Summer 2026). I'm using this to keep track of my code as I work through the project — each commit roughly corresponds to a step or a meaningful chunk of progress.

## What the project is about

A fire authority manages a 21×21 grid landscape and needs to decide where to place firebreaks (cleared cells that block fire). There's a settlement in the NE corner and a wetland in the NW corner that need protecting. The total budget is 17 cells: 13 placed before the fire season, and 4 placed reactively after wind conditions are revealed.

The project walks through 4 steps of increasing complexity:

1. **Step 1** — Model the grid as a network and use max-flow to find where fire pressure concentrates (bottleneck analysis)
2. **Step 2** — Formulate budget-constrained MILPs (two models: binary propagation and continuous reachability) to optimise firebreak placement under deterministic conditions
3. **Step 3** — Evaluate the Step 2 plan under wind/ignition uncertainty using Monte Carlo simulation. Need to build an acceptance-rejection sampler for the wind direction (mixed von Mises)
4. **Step 4** — Introduce reactive firebreaks: solve a per-scenario MILP after wind is revealed, compare against purely permanent plans

## My approach

I'm planning to do Steps 1–3 thoroughly and give Step 4 a solid attempt. The main tools are AMPL (for optimisation models) and R (for simulation and data processing). I call AMPL from R via `system2()`.

For the report, I'm writing directly on Cadmus. This repo is just for the code side of things.

## Environment

- R 4.5.1 on macOS (Apple Silicon)
- AMPL v20250901 with academic license
- Solvers: Gurobi 13.0 (primary), HiGHS 1.11, CBC (backups)
- `.Rprofile` sets up all the solver paths automatically — just `setwd()` to this directory and everything should work

## File structure (planned)

```
.Rprofile              <- auto-configures AMPL/solver paths
build_dat.R            <- generates .dat files from landscape data
step1_maxflow.mod/dat/run
step2_model_a.mod/dat/run
step2_model_b.mod/dat/run
step3_mc.R             <- AR sampler, MC evaluation, risk measures
step4_reactive.mod
step4_ramp.R           <- per-scenario reactive solve loop
figures/               <- output plots
scratch.md             <- my working notes (not submitted)
```

## Note on AI usage

I used Claude (education version) to help draft some initial code, which I then adapted and debugged to fit my specific data and model requirements. The full details are in my AI declaration in the submitted report.
