# Model B key differences from Model A:

param alpha  default 0.5;              # compression exponent (< 1)

# Derived parameter: compressed edge weight
param w {(u,v) in EDGES} = p_burn[u,v] ^ alpha;

# Continuous fire reachability instead of binary burn
var x     {CELLS} >= 0, <= 1;

# Objective: minimise total weighted reachability
minimize Total_Damage:
    sum {t in TARGETS} weight[t] * x[t];

# Ignition: row-21 cells have full reachability
subject to Ignition_Reach {u in IGNITION}:
    x[u] = 1;

# A cleared cell has zero reachability
subject to Clear_Blocks {v in CELLS}:
    x[v] <= 1 - clear[v];

# Continuous propagation: reachability flows multiplicatively
#    x[v] >= w[u,v] * x[u] - clear[v]
#
#    If v not cleared (clear[v]=0): x[v] >= w * x[u]
#    If v cleared     (clear[v]=1): x[v] >= w * x[u] - 1
#                                   which is non-binding since x >= 0
subject to Propagation {(u, v) in EDGES}:
    x[v] >= w[u,v] * x[u] - clear[v];
