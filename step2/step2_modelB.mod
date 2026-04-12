# ============================================================
# Step 2 — Model B: Continuous reachability firebreak MILP
# ============================================================
# Minimise weighted damage to targets by clearing at most B
# cells. Fire reachability propagates multiplicatively along
# edges with compressed weights w = p_burn^alpha.
#
# At optimality x[v] equals the product of compressed edge
# weights along the strongest unblocked path from the ignition
# front to v. The compression (alpha < 1) prevents values from
# decaying below solver tolerance over many hops.
# ============================================================

# --- Sets ---------------------------------------------------
set CELLS;                              # all vegetated cells ("r_c")
set IGNITION within CELLS;              # row-21 cells that start on fire
set TARGETS  within CELLS;              # settlement + wetland cells
set EDGES    within {CELLS, CELLS};     # directed spread edges

# --- Parameters ---------------------------------------------
param weight {TARGETS};                 # damage weight per target
param p_burn {EDGES};                   # raw spread probability
param B      >= 0, default 13;          # clearing budget
param alpha  default 0.5;              # compression exponent (< 1)

# --- Derived parameter: compressed edge weight ---------------
param w {(u,v) in EDGES} = p_burn[u,v] ^ alpha;

# --- Decision variables -------------------------------------
var clear {CELLS} binary;               # 1 if cell is cleared
var x     {CELLS} >= 0, <= 1;           # continuous fire reachability

# --- Objective: minimise total weighted reachability ---------
minimize Total_Damage:
    sum {t in TARGETS} weight[t] * x[t];

# --- Constraints --------------------------------------------

# 1. Budget: clear at most B cells
subject to Budget:
    sum {v in CELLS} clear[v] <= B;

# 2. Ignition: row-21 cells have full reachability, never cleared
subject to Ignition_Reach {u in IGNITION}:
    x[u] = 1;

subject to Ignition_NoClear {u in IGNITION}:
    clear[u] = 0;

# 3. Targets: cannot be cleared
subject to Target_NoClear {t in TARGETS}:
    clear[t] = 0;

# 4. A cleared cell has zero reachability
subject to Clear_Blocks {v in CELLS}:
    x[v] <= 1 - clear[v];

# 5. Continuous propagation: reachability flows multiplicatively
#    x[v] >= w[u,v] * x[u] - clear[v]
#
#    If v not cleared (clear[v]=0): x[v] >= w * x[u]
#    If v cleared     (clear[v]=1): x[v] >= w * x[u] - 1
#                                   which is non-binding since x >= 0
subject to Propagation {(u, v) in EDGES}:
    x[v] >= w[u,v] * x[u] - clear[v];
