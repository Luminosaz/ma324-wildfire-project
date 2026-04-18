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
