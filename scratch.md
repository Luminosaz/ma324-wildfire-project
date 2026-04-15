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
- Pilot N=50, alpha=0.1, R=4, Gurobi timelim=120s mipgap=0.01
  - none:      E=13.50, CVaR90=74.00
  - perm_13:   E=3.50,  CVaR90=21.00
  - perm_17:   E=2.60,  CVaR90=13.75
  - reactive:  E=8.90,  CVaR90=39.38   ← 反常地高
- Runtime N=50 ≈ 636s (~12s/solve). N=1000 ≈ 3.5 h.
- Per-scenario diff (reactive − perm_13): mean +5.40, **median 0**
  → 大多数 scenario reactive 不输,少数被 RNG 拉高

### Discussion bullets
- (fill in)

### 待问老师的问题（Step 4 邮件合集）

1. **CRN 的 RNG 策略**（最关键,影响 Q2 的整个分析）

   **背景**: Step 3 `run_mc_evaluation` 用全局 `set.seed(1)` 一次,之后 for plan × for scenario 顺序调 `simulate_fire`。simulator 内部的 `runif` 在 plan 之间自然推进,4 个 plan **不共享 simulator random stream**,只共享 (ignition, V, θ)。

   **Step 4 的问题**: reactive 的 firebreak set = `perm_13 ∪ react_cells`(严格 superset)。物理上在**同一 simulator realisation** 下必然 `damage(reactive) ≤ damage(perm_13)`。但 Step 3 式 seed 做法导致单 scenario 层面经常出现 `damage(reactive) > damage(perm_13)`——这只可能来自两个 plan 的 simulator runif 序列不同。

   Q2 的核心是 per-scenario "reactive effort gap",这个 gap 在当前 seed 做法下被 RNG 噪声污染(mean +5.40 但 median 0)。

   **想做的修法**(per-scenario derived seed, 仍然从全局 set.seed(1) 单一派生):
   ```r
   for i in 1..N:
     sim_seed = 1000L + i
     set.seed(sim_seed); simulate_fire(land_none, ...)
     set.seed(sim_seed); simulate_fire(land_13,   ...)
     set.seed(sim_seed); simulate_fire(land_17,   ...)
     # solve_reactive 不动 RNG
     set.seed(sim_seed); simulate_fire(land_react, ...)
   ```

   **要老师确认的两问**:
   (a) Step 4 是否允许采用 per-scenario derived seed(`set.seed(1000 + i)`),让 4 个 plan 共享 scenario-level simulator stream?
   (b) 若不允许,Step 3 现有做法导致 Q2 的 per-scenario gap 分析被 RNG 污染,这在 Discussion 里怎么 frame 合适?
   (c) Step 3 的数值结果是否需要回头用更严格的 CRN 重跑(影响 proposed vs corridor 的 4.29 vs 3.89 比较)?

2. **B=17 的 MIP gap 问题**（Q3 素材,已部分自答,但想确认解读）
   - Model B B=17 α=0.1 在 120s 下 reported relmipgap=0.254(25%),1800s 下达 optimal,objective 都是 11.585。
   - 解读:120s 的 incumbent 其实就是 global optimal,gap 25% 仅是 Gurobi 的 LP 下界没 tighten,不代表解次优。
   - 要确认:Discussion Q3 里把"gap ≠ solution quality"这个解读写出来,合理吗?还是应该更保守地说"gap 大意味着我们的 reactive 评估可能 tail 悲观"?

3. **Ignition proxy 在 reactive MILP 里的定义**（已按 brief 推荐用 southern-front row 21,想确认）
   - 目前 reactive MILP 的 IGNITION = row 21 所有 vegetated cells(Step 2 的做法沿用)。
   - brief 原文 Q6 也明确说"reactive MILP fixes ignition as a southern-front proxy"。
   - 但真实 ignition 是 uniform over vegetated cells(N=50 下有 ~10 个 scenario ignition 在 row 10+)。这个 mismatch 正是 Q2 的 dominance violation 的主要嫌疑。
   - 要确认:southern-front proxy 是 Step 4 **唯一允许的 planning-stage ignition assumption**,还是可以讨论 alternative(如 expected-over-ignition-distribution)作为 depth extension?
