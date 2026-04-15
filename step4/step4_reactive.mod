# ============================================================
# Step 4 — Reactive firebreak MILP (Model B extended)
# ============================================================
# Given a fixed permanent plan (set PERM from Step 3) and a
# revealed wind scenario, choose at most R additional reactive
# cells to minimise weighted target reachability.
#
# Differences from Step 2 Model B:
#   - No `var clear` or `param B`.  Permanent cells are DATA
#     (set PERM), reactive cells are the sole decision variable.
#   - Edge weights `w` are pre-computed in R by edges_wind()
#     for the revealed wind, replacing p_burn^alpha.
#   - A cell is blocked iff it is permanent OR reactive.
# ============================================================

# --- Sets ---------------------------------------------------
set CELLS;                              # all vegetated cells ("r_c")
set IGNITION within CELLS;              # row-21 cells that start on fire
set TARGETS  within CELLS;              # settlement + wetland cells
set EDGES    within {CELLS, CELLS};     # directed spread edges
set PERM     within CELLS;              # permanently cleared cells (fixed)

# --- Data-integrity checks ----------------------------------
check {u in IGNITION}: u not in PERM;
check {t in TARGETS}:  t not in PERM;

# --- Parameters ---------------------------------------------
param weight {TARGETS};                 # damage weight per target
param w      {EDGES};                   # pre-computed edge weight (p_burn^alpha)
param R      >= 0;                      # reactive clearing budget

# --- Decision variables -------------------------------------
var react {CELLS} binary;               # 1 if cell is reactively cleared
var x     {CELLS} >= 0, <= 1;           # continuous fire reachability

# --- Objective: minimise total weighted reachability ---------
minimize Total_Damage:
    sum {t in TARGETS} weight[t] * x[t];

# --- Constraints --------------------------------------------

# 1. Reactive budget: clear at most R additional cells
subject to Reactive_Budget:
    sum {v in CELLS} react[v] <= R;

# 2. Ignition: full reachability, never cleared
subject to Ignition_Reach {u in IGNITION}:
    x[u] = 1;

subject to Ignition_NoReact {u in IGNITION}:
    react[u] = 0;

# 3. Targets: cannot be reactively cleared
subject to Target_NoReact {t in TARGETS}:
    react[t] = 0;

# 4. Permanent cells: already cleared — zero reachability,
#    no reactive spend
subject to Perm_Reach {v in PERM}:
    x[v] = 0;

subject to Perm_NoReact {v in PERM}:
    react[v] = 0;

# 5. A reactively cleared cell has zero reachability
subject to React_Blocks {v in CELLS diff PERM}:
    x[v] <= 1 - react[v];

# 6. Continuous propagation with wind-dependent weights
#    x[v] >= w[u,v] * x[u] - (blocked_v)
#
#    For v in PERM:  x[v] = 0 already (constraint 4), so
#                    propagation is non-binding.
#    For u in PERM:  x[u] = 0, so RHS = w * 0 - react[v] <= 0,
#                    non-binding since x[v] >= 0.
#    For u, v not in PERM:
#      react[v] = 0 → x[v] >= w * x[u]
#      react[v] = 1 → x[v] >= w * x[u] - 1  (non-binding)
subject to Propagation {(u, v) in EDGES: u not in PERM and v not in PERM}:
    x[v] >= w[u,v] * x[u] - react[v];
