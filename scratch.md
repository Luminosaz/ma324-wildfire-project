# MA324 Project Lab Notebook (scratch — not submitted)

> Daily numbers, parameter choices, and discussion bullets. Used for Cadmus writing reference. **Not for submission.**

---

## Environment Setup (2026-04-08)

**AMPL binary** (via amplpy): `/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_base/bin/ampl` (v20250901)
**Gurobi solver**: `/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_gurobi/bin/gurobi` (v13.0.0)
**HiGHS solver**: `/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_highs/bin/highs` (v1.11.0)
**CBC solver**: `/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_cbc/bin/cbc`
**License**: AMPL Academic (j.zhang227@lse.ac.uk), no size limits, expires 2026-05-10
**License file**: `/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_base/bin/ampl.lic`
**R version**: 4.5.1 (aarch64-apple-darwin20)

### Setup verified 2026-04-08
- AMPL binary: ✅
- HiGHS solver: ✅ (tested with 3585-constraint MILP, optimal)
- Gurobi solver: ✅ (tested, optimal)
- R → AMPL via system2(): ✅
- .Rprofile sets AMPL_BIN / GUROBI_BIN / HIGHS_BIN / CBC_BIN and adds to PATH

---

## Day 1 (2026-04-10) — Step 1 max-flow

### Numbers
- Settlement max flow = 4.06
- Wetland max flow = 3.795
- Total edges (no wind): 3248
- Vegetated cells: 439, bare ground: 2

### Discussion bullets
- Settlement 比 wetland 更难保护（max flow 4.06 vs 3.795），火压力更大
- 瓶颈集中在目标区域紧邻外围（settlement: cols 15-16 + row 7, wetland: cols 7-8 + row 7），不在 current corridor (row 10)
- Current corridor 放太远了——row 10 离两个目标都有 4 行距离，火可以绕过去。bottleneck 分析说防线应该贴着目标放

---

## Day 2 (2026-04-11/12) — Step 2 Models

### Numbers
- Model A (τ=0, B=13): damage=80, clears col15(r1-7)+row7(c16-21)=13 cells, protects settlement only
- Model A: τ≤0.1 → all targets burn (240), τ≥0.105 → no fire reaches targets (0). No gradual transition.
- Model B (α=0.1, B=13): damage=14.07, same col15+row7 pattern
- Model B (α=0.1, B=16): damage=11.98, adds wetland protection (11_4, 11_5, 11_6)
- Model B (α=0.05, B=5): damage=95.1 → sim=10.1
- Model B (α=0.05, B=8): damage=87.4 → sim=4.8
- Baseline (no firebreaks): sim mean=16.4 (K=1000)
- Current corridor (row10, 17 cells): sim mean=1.6 (K=1000)
- Model A B=13 τ=0: sim mean=4.8 (K=1000)
- Model B B=16 α=0.05: sim mean=3.0 (K=1000)
- B=20 (both models): sim=0, full protection

### Discussion bullets
- Model A 是 all-or-nothing: τ 控制的是"火能不能到 target"的开关，没有中间地带。不适合做 budget sensitivity
- Model B 有平滑的 budget-damage 曲线，能回答"多一格值不值"
- 两个模型都严重高估 damage（Model A: 80 vs sim 4.8, Model B α=0.05: 95 vs sim 10.1），但 plan 质量不错——model can predict poorly but produce a good plan
- B=13 时两个模型选了相同的 firebreak 位置（col15+row7），都只保护 settlement
- B=16 时 Model B 开始分配资源给 wetland（加了 11_4/5/6），Model A 仍然只保护 settlement（只用 13 格子）
- Current corridor (17 cells) sim=1.6 其实比 MILP B=13 (sim=4.8) 还好，但 corridor 用了更多格子
- α=0.5 太大完全不能用（reachability 15 hop 后趋近 0），α=0.1 左右合适
- marginal value: Model B 在 B=10→13 有最大跳跃（settlement 被保护住的临界点）

### 待问老师的问题（Step 2 邮件合集）
1. **Model A τ 的选择**: binary propagation 在 τ=0.10 和 τ=0.105 之间有 cliff（240→0），没有中间值。这是预期的结果还是我的 formulation 有问题？
2. **Model B α 的选择**: α=0.5 时 reachability 在 15+ hops 后衰减到 <0.001，模型失去区分度。α=0.1 效果好但 brief 原文说 "α < 1"，没说具体范围。用 α=0.1 这么小的值合理吗？
3. **Model A budget 不生效**: τ=0 时 B=5/8/10 的最优解是不清任何格子（damage=240），因为小预算无法切断所有路径。模型返回 n_cleared=0 是正确行为还是 bug？
4. **Shadow price on budget**: brief 说 "obtain the shadow price on the budget constraint"，但 MILP 没有 LP dual。我用 marginal value = damage(B) - damage(B+1) 来近似，这个做法可以吗？
5. **代码提交格式**: 报告里贴完整 R 脚本（~200 行 step2_simulate.R）会很占字数。可以只贴核心函数+调用骨架，还是必须完整？

---

## Day 3 (2026-04-11) — Step 3 MC

### Numbers
(fill in)

### Discussion bullets
- (fill in)

---

## Day 4 (2026-04-12) — Step 4 reactive

### Numbers
(fill in)

### Discussion bullets
- (fill in)
