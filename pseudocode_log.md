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

(todo)

### model A: binary propagation

(todo)

### model B: continuous reachability

(todo)

---

## step 3 — monte carlo

### AR sampler for wind direction

(todo)

### MC evaluation

(todo)

---

## step 4 — reactive

(todo)
