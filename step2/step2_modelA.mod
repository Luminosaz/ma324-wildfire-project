# ============================================================
# Step 2 — Model A: Binary propagation firebreak MILP
# ============================================================
# Minimise weighted damage to targets by clearing at most B
# cells. Fire propagates deterministically along edges whose
# p_burn exceeds threshold τ (filtered inside the model).
# ============================================================

# --- Sets ---------------------------------------------------
set CELLS;                              # all vegetated cells ("r_c")
set IGNITION within CELLS;              # row-21 cells that start on fire
set TARGETS  within CELLS;              # settlement + wetland cells
set EDGES    within {CELLS, CELLS};     # all directed spread edges

# --- Parameters ---------------------------------------------
param weight {TARGETS};                 # damage weight per target
param p_burn {EDGES};                   # spread probability
param B      >= 0, default 13;          # clearing budget
param tau    default 0;                 # edge activation threshold

# --- Decision variables -------------------------------------
var clear {CELLS} binary;               # 1 if cell is cleared
var burn  {CELLS} binary;               # 1 if cell burns

# --- Objective: minimise total weighted damage ---------------
minimize Total_Damage:
    sum {t in TARGETS} weight[t] * burn[t];

# --- Constraints --------------------------------------------

# 1. Budget: clear at most B cells
subject to Budget:
    sum {v in CELLS} clear[v] <= B;

# 2. Ignition: row-21 cells always burn, never cleared
subject to Ignition_Burns {u in IGNITION}:
    burn[u] = 1;

subject to Ignition_NoClear {u in IGNITION}:
    clear[u] = 0;

# 3. Targets: cannot be cleared (they must be protected by blocking fire earlier)
subject to Target_NoClear {t in TARGETS}:
    clear[t] = 0;

# 4. A cleared cell cannot burn
subject to Clear_Blocks {v in CELLS}:
    burn[v] <= 1 - clear[v];

# 5. Binary propagation: if u burns and v is not cleared, v burns
#    Only active edges (p_burn > τ) generate constraints.
#    burn[v] >= burn[u] - clear[v]   for each active edge (u,v)
subject to Propagation {(u, v) in EDGES: p_burn[u,v] > tau}:
    burn[v] >= burn[u] - clear[v];
