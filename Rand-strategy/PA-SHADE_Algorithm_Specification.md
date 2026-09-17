# PA-SHADE: Phase-Adaptive Success-History based Adaptive Differential Evolution

## Algorithm Specification and Pseudo-code

---

## 1. 算法概述

**PA-SHADE** (Phase-Adaptive SHADE) 是在经典 SHADE (Success-History based Adaptive Differential Evolution) 框架基础上，引入**三阶段动态自适应策略**的改进算法：

- **探索期** (FEs < 0.6·maxFEs): 大范围搜索，保持种群多样性
- **过渡期** (0.25·maxFEs ≤ FEs < 0.6·maxFEs): 线性混合交叉，平衡探索与利用  
- **收敛期** (FEs ≥ 0.6·maxFEs): 适应度导向变异，定向收敛

---

## 2. 核心符号表

| 符号 | 含义 | 默认值 |
|:---:|:---|:---|
| $NP$ | 初始种群规模 | 1018 (D=50) |
| $D$ | 问题维度 | 50 |
| $maxFEs$ | 最大函数评估次数 | $10000 \times D$ |
| $H$ | 历史记忆库大小 | 5 |
| $A_{init}$ | 初始存档容量上限 | $NP \times 1.5$ |
| $pb$ | pbest 比例 | 0.085 ~ 0.17 (动态) |
| $\mathcal{F}_i$ | 第 $i$ 个个体的缩放因子 | 自适应 |
| $Cr_i$ | 第 $i$ 个个体的交叉率 | 自适应 |
| $M_\mathcal{F}$ | 历史记忆中 $\mathcal{F}$ 的均值 | 自适应更新 |
| $M_{Cr}$ | 历史记忆中 $Cr$ 的均值 | 自适应更新 |
| $dy\_NP$ | 当前种群规模 (动态递减) | $NP \rightarrow NP_{min}$ |
| $archive$ | 外部存档集合 | 动态维护 |

---

## 3. 算法流程图

```
┌─────────────────────────────────────────────┐
│              初始化阶段                        │
│  • 随机初始化种群 X[1..NP][1..D]             │
│  • 评估适应度 f[1..NP]                      │
│  • 初始化历史记忆 his_mem[1..H] = (0.3, 0.8) │
│  • 外部存档 archive ← ∅                       │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│         主循环 (while FEs < maxFEs)          │
└─────────────────────────────────────────────┘
                    ↓
        ┌───────────────────┐
        │   1. pbest 排序    │
        │  按适应度升序排列   │
        │  pbSize = max(2, ⌊pb·dy_NP⌋) │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   2. 参数自适应    │
        │  • 读取历史记忆 (MF, Mcr) │
        │  • 柯西分布采样 F  │
        │  • 高斯分布采样 Cr │
        │  • 分阶段调整 F 下限│
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   3. 变异操作      │
        │  • 选择 pbest 个体 │
        │  • 动态邻域选择 r2 │
        │  • 存档扩展选择 r3 │
        │  • 分阶段变异策略  │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   4. 交叉操作      │
        │  • 二项式交叉      │
        │  • 线性混合:      │
        │    U = (1-c)·V + c·X │
        │    c = 0.3·FEs/maxFEs │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   5. 选择操作      │
        │  • 贪婪选择        │
        │  • 更新外部存档    │
        │  • 记录成功参数    │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   6. 历史记忆更新  │
        │  • 加权 Lehmer 平均 F │
        │  • 加权算术平均 Cr │
        │  • 环形缓冲区索引  │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   7. 种群规模调整  │
        │  • 线性递减 dy_NP  │
        │  • 移除最差个体    │
        │  • 调整存档上限    │
        └───────────────────┘
                    ↓
        ┌───────────────────┐
        │   8. pbest 比例调整│
        │  pb = 0.085 + 0.085·FEs/maxFEs │
        └───────────────────┘
                    ↓
              [循环结束?]
                    ↓
        ┌───────────────────┐
        │      输出最优解    │
        └───────────────────┘
```

---

## 4. 核心创新组件详解

### 4.1 分阶段变异策略 (Phase-Dependent Mutation)

| 阶段 | 条件 | 变异公式 | 目的 |
|:---:|:---|:---|:---|
| **探索期** | $FEs < 0.6 \cdot maxFEs$ | $V_i = X_{r_4} + \mathcal{F}_i(X_{pbest} - X_{r_4}) + \mathcal{F}_i(X_{r_2} - X_{r_3})$ | 广泛搜索，保持多样性 |
| **收敛期** | $FEs \geq 0.6 \cdot maxFEs$ | $V_i = X_{r_4} + \mathcal{F}_i(X_{pbest} - X_{r_4}) + \mathcal{F}_i(X_{better} - X_{worse})$ | 定向收敛，利用适应度信息 |

其中 $X_{better}$ 和 $X_{worse}$ 根据 $f(r_2)$ 和 $f(r_3)$ 动态确定：

$$
\text{if } f(r_2) < f(r_3): \quad better = r_2, \; worse = r_3 \\text{else}: \quad better = r_3, \; worse = r_2
$$

### 4.2 动态邻域范围选择 (Dynamic Neighborhood Range)

| 阶段 | neighbour₂ 选择范围 | 目的 |
|:---:|:---|:---|
| 早期 | $np\_range = \max(3, \lfloor 0.7 \times dy\_NP \rfloor)$ | 大范围优质个体采样 |
| 后期 | $np\_range = \max(3, \lfloor (0.5 - 0.2 \frac{FEs}{maxFEs}) \times dy\_NP \rfloor)$ | 逐步收缩至前 30% |

### 4.3 线性混合交叉 (Linear Blending Crossover)

混合系数随进化进程线性增长：

$$
c = 0.3 \times \frac{FEs}{maxFEs}
$$

试验向量生成：

$$
U_{ij} = (1-c) \cdot V_{ij} + c \cdot X_{ij}
$$

- **早期** ($c \approx 0$): 主要使用变异向量 $V$，强调探索
- **后期** ($c \approx 0.3$): 逐渐保留原个体信息 $X$，强调利用

### 4.4 存档替换策略 (Archive Replacement with Worst Pool)

存档满时，从**最差 30% 池中随机选择**替换：

1. 对存档按适应度排序，找出最差 30% 个体
2. 在该池内随机选择替换目标
3. 避免总是替换全局最差，维持存档多样性

---

## 5. 算法伪代码

```
Algorithm: PA-SHADE (Phase-Adaptive SHADE)
Input:  Population size NP, Dimension D, Max FEs maxFEs, 
        History memory size H, Archive limit factor
Output: Best solution X_best, Best fitness f_best

─────────────────────────────────────────────
// 1. INITIALIZATION
─────────────────────────────────────────────
for i = 1 to NP do
    for j = 1 to D do
        X[i][j] ← Random(lb, ub)
    end for
    f[i] ← Evaluate(X[i])
    FEs ← FEs + 1
end for

for k = 1 to H do
    his_mem[k] ← (0.3, 0.8)       // (MF, Mcr)
end for

archive ← ∅
his_f ← 0                          // Current history index
his_k ← 1                          // Next update index
pb ← 0.085
dy_NP ← NP
ANum ← ⌊NP × 1.5⌋

─────────────────────────────────────────────
// 2. MAIN LOOP
─────────────────────────────────────────────
while FEs < maxFEs do

    // ── 2.1 PBEST SORTING ──
    bestp ← SortIndicesByFitness({1, 2, ..., dy_NP})
    pbSize ← max(2, ⌊pb × dy_NP⌋)

    // ── 2.2 PARAMETER ADAPTATION ──
    for i = 1 to dy_NP do
        MF ← his_mem[his_f].MF
        Mcr ← his_mem[his_f].Mcr

        // Sample scaling factor F using Cauchy distribution
        F_i ← Cauchy(MF, 0.1)
        while F_i ≤ 0 do
            F_i ← Cauchy(MF, 0.1)
        end while

        // Phase-dependent F floor adjustment
        if FEs < 0.25 × maxFEs and F_i < 0.5 then
            F_i ← 0.5
        end if
        if F_i > 1 then F_i ← 1 end if

        // Sample crossover rate Cr using Gaussian distribution
        if Mcr = -1 then
            Cr_i ← 0
        else
            Cr_i ← Gaussian(Mcr, 0.1)
            if FEs < 0.25 × maxFEs then
                Cr_i ← max(Cr_i, 0.5)
            end if
        end if
        Cr_i ← Clip(Cr_i, 0, 1)
    end for

    // ── 2.3 MUTATION ──
    for i = 1 to dy_NP do
        // Select pbest individual
        r₁ ← Random(1, pbSize)
        pbest ← bestp[r₁]
        while pbest = i do
            r₁ ← Random(1, pbSize)
            pbest ← bestp[r₁]
        end while

        // Dynamic neighborhood range for r₂
        if FEs < 0.5 × maxFEs then
            np_range ← max(3, ⌊0.7 × dy_NP⌋)
        else
            np_range ← max(3, ⌊(0.5 - 0.2 × FEs/maxFEs) × dy_NP⌋)
        end if

        r₂ ← Random(1, np_range)
        nr₂ ← bestp[r₂]
        while nr₂ = i or nr₂ = pbest do
            r₂ ← Random(1, np_range)
            nr₂ ← bestp[r₂]
        end while

        // Select r₃ from population + archive
        r₃_range ← dy_NP + |archive|
        r₃ ← Random(1, r₃_range)
        while r₃ = i or r₃ = pbest or r₃ = nr₂ do
            r₃ ← Random(1, r₃_range)
        end while
        if r₃ > dy_NP then
            nr₃ ← archive[r₃ - dy_NP]
        else
            nr₃ ← r₃
        end if

        // Select random base vector r₄
        r₄ ← Random(1, dy_NP)
        while r₄ ∈ {i, pbest, nr₂, nr₃} do
            r₄ ← Random(1, dy_NP)
        end while

        // Phase-dependent mutation strategy
        if FEs < 0.6 × maxFEs then
            // Exploration phase: standard DE/rand/1 variant
            for j = 1 to D do
                V[i][j] ← X[r₄][j] + F_i × (X[pbest][j] - X[r₄][j])
                            + F_i × (X[nr₂][j] - X[nr₃][j])
            end for
        else
            // Convergence phase: fitness-directed mutation
            if f[nr₂] < f[nr₃] then
                better ← nr₂;  worse ← nr₃
            else
                better ← nr₃;  worse ← nr₂
            end if
            for j = 1 to D do
                V[i][j] ← X[r₄][j] + F_i × (X[pbest][j] - X[r₄][j])
                            + F_i × (X[better][j] - X[worse][j])
            end for
        end if

        // Boundary handling
        for j = 1 to D do
            while V[i][j] < lb do
                V[i][j] ← (lb + X[i][j]) / 2
            end while
            while V[i][j] > ub do
                V[i][j] ← (ub + X[i][j]) / 2
            end while
        end for
    end for

    // ── 2.4 CROSSOVER (Linear Blending) ──
    c ← 0.3 × (FEs / maxFEs)              // Blending coefficient

    for i = 1 to dy_NP do
        j_rand ← Random(1, D)              // Mandatory inheritance dimension
        for j = 1 to D do
            if Random(0,1) ≤ Cr_i or j = j_rand then
                U[i][j] ← (1-c) × V[i][j] + c × X[i][j]    // Linear blending
            else
                U[i][j] ← X[i][j]
            end if
        end for
    end for

    // ── 2.5 SELECTION + ARCHIVE UPDATE ──
    S_F ← ∅,  S_Cr ← ∅,  S_Δf ← ∅        // Successful parameter sets

    for i = 1 to dy_NP and FEs < maxFEs do
        fU ← Evaluate(U[i])
        FEs ← FEs + 1

        if fU < f[i] then                  // Strict improvement
            // Store in archive
            if |archive| < ANum then
                archive ← archive ∪ {(X[i], f[i])}
            else
                // Replace from worst 30% pool
                worst_pool ← Bottom30%(archive)
                idx ← Random(worst_pool)
                archive[idx] ← (X[i], f[i])
            end if

            S_F ← S_F ∪ {F_i}
            S_Cr ← S_Cr ∪ {Cr_i}
            S_Δf ← S_Δf ∪ {f[i] - fU}
        end if

        // Greedy selection (≤ allows equal replacement)
        if fU ≤ f[i] then
            X[i] ← U[i]
            f[i] ← fU
        end if
    end for

    // ── 2.6 HISTORY MEMORY UPDATE ──
    if |S_F| > 0 then
        sum_Δf ← Σ_{k=1}^{|S_F|} S_Δf[k]

        // Weighted mean for Cr
        t₁ ← Σ (S_Δf[k] / sum_Δf) × S_Cr[k]
        t₂ ← Σ (S_Δf[k] / sum_Δf) × S_Cr[k]²
        if t₁ ≠ 0 then
            mean_Cr ← t₂ / t₁
            his_mem[his_k].Mcr ← (mean_Cr + his_mem[his_f].Mcr) / 2
        else
            his_mem[his_k].Mcr ← -1
        end if

        // Weighted Lehmer mean for F
        k₁ ← Σ (S_Δf[k] / sum_Δf) × S_F[k]²
        k₂ ← Σ (S_Δf[k] / sum_Δf) × S_F[k]
        mean_F ← k₁ / k₂
        his_mem[his_k].MF ← (mean_F + his_mem[his_f].MF) / 2

        his_k ← (his_k + 1) mod H
        his_f ← (his_f + 1) mod H
    end if

    // ── 2.7 POPULATION SIZE REDUCTION ──
    new_NP ← round(NP + (NP_min - NP) × FEs / maxFEs)
    if new_NP < dy_NP then
        while dy_NP > new_NP do
            worst_idx ← argmax{f[i] | i ∈ [1, dy_NP]}
            Swap(X[worst_idx], X[dy_NP])
            Swap(f[worst_idx], f[dy_NP])
            dy_NP ← dy_NP - 1
        end while
    end if

    // Adjust archive limit
    ANum_new ← ⌊dy_NP × 1.5⌋
    if |archive| > ANum_new then
        while |archive| > ANum_new do
            worst_a ← argmax{f[a] | a ∈ archive}
            Remove(archive, worst_a)
        end while
    end if
    ANum ← ANum_new

    // ── 2.8 PBEST RATIO ADJUSTMENT ──
    pb ← 0.085 + 0.085 × (FEs / maxFEs)

end while

return X[argmin_i f[i]], min_i f[i]
```

---

## 6. 关键公式速查表

### 6.1 参数采样

| 参数 | 分布 | 公式 | 截断处理 |
|:---:|:---:|:---|:---|
| $\mathcal{F}_i$ | 柯西分布 | $\mathcal{F}_i \sim \text{Cauchy}(M_\mathcal{F}, 0.1)$ | 重采样至 $\mathcal{F}_i > 0$，上限为 1 |
| $Cr_i$ | 高斯分布 | $Cr_i \sim \mathcal{N}(M_{Cr}, 0.1)$ | 截断至 $[0, 1]$ |

### 6.2 分阶段变异

**探索期** ($FEs < 0.6 \cdot maxFEs$):

$$
\mathbf{V}_i = \mathbf{X}_{r_4} + \mathcal{F}_i \cdot (\mathbf{X}_{pbest} - \mathbf{X}_{r_4}) + \mathcal{F}_i \cdot (\mathbf{X}_{r_2} - \mathbf{X}_{r_3})
$$

**收敛期** ($FEs \geq 0.6 \cdot maxFEs$):

$$
\mathbf{V}_i = \mathbf{X}_{r_4} + \mathcal{F}_i \cdot (\mathbf{X}_{pbest} - \mathbf{X}_{r_4}) + \mathcal{F}_i \cdot (\mathbf{X}_{better} - \mathbf{X}_{worse})
$$

其中：

$$
(better, worse) = \begin{cases} (r_2, r_3) & \text{if } f(r_2) < f(r_3) \\ (r_3, r_2) & \text{otherwise} \end{cases}
$$

### 6.3 线性混合交叉

$$
c = 0.3 \times \frac{FEs}{maxFEs}
$$

$$
U_{ij} = \begin{cases} (1-c) \cdot V_{ij} + c \cdot X_{ij} & \text{if } rand \leq Cr_i \text{ or } j = j_{rand} \\ X_{ij} & \text{otherwise} \end{cases}
$$

### 6.4 历史记忆更新 (加权 Lehmer 平均)

权重：

$$
w_k = \frac{\Delta f_k}{\sum_{m} \Delta f_m}, \quad \Delta f_k = f_{old} - f_{new}
$$

$\mathcal{F}$ 更新：

$$
\text{mean}_\mathcal{F} = \frac{\sum w_k \mathcal{F}_k^2}{\sum w_k \mathcal{F}_k}, \quad M_\mathcal{F}^{new} = \frac{\text{mean}_\mathcal{F} + M_\mathcal{F}^{old}}{2}
$$

$Cr$ 更新：

$$
\text{mean}_{Cr} = \frac{\sum w_k Cr_k^2}{\sum w_k Cr_k}, \quad M_{Cr}^{new} = \frac{\text{mean}_{Cr} + M_{Cr}^{old}}{2}
$$

### 6.5 动态调整

| 组件 | 公式 | 范围 |
|:---:|:---|:---|
| 种群规模 | $NP_{new} = NP + (NP_{min} - NP) \times \frac{FEs}{maxFEs}$ | $NP \rightarrow NP_{min}$ |
| pbest 比例 | $pb = 0.085 + 0.085 \times \frac{FEs}{maxFEs}$ | $0.085 \rightarrow 0.17$ |
| 混合系数 | $c = 0.3 \times \frac{FEs}{maxFEs}$ | $0 \rightarrow 0.3$ |
| 邻域范围 | $np\_range = \max(3, \lfloor (0.5 - 0.2 \frac{FEs}{maxFEs}) \times dy\_NP \rfloor)$ | $0.7\cdot NP \rightarrow 0.3\cdot NP$ |

---

## 7. 与标准 SHADE 的对比

| 特性 | SHADE | PA-SHADE (本算法) |
|:---:|:---|:---|
| **变异策略** | 固定 DE/current-to-pbest/1 | **分阶段**：前期标准变异 / 后期适应度导向 |
| **交叉方式** | 二项式交叉 | **线性混合交叉**，动态平衡探索与利用 |
| **邻域范围** | 固定 pbest 比例 | **动态收缩**，从 70% 逐步降至 30% |
| **存档替换** | 随机替换 | **最差 30% 池随机替换**，维持多样性 |
| **pbest 比例** | 固定值 | **线性递增**，从 0.085 增至 0.17 |
| **F 下限** | 固定 | **分阶段**：前期强制 $\mathcal{F} \geq 0.5$ |

---

## 8. 复杂度分析

| 操作 | 时间复杂度 | 空间复杂度 |
|:---|:---:|:---:|
| 初始化 | $O(NP \times D)$ | $O(NP \times D)$ |
| pbest 排序 | $O(dy\_NP \log dy\_NP)$ | $O(dy\_NP)$ |
| 参数自适应 | $O(dy\_NP)$ | $O(1)$ |
| 变异操作 | $O(dy\_NP \times D)$ | $O(dy\_NP \times D)$ |
| 交叉操作 | $O(dy\_NP \times D)$ | $O(dy\_NP \times D)$ |
| 选择 + 存档 | $O(dy\_NP)$ | $O(ANum)$ |
| 历史记忆更新 | $O(|S|)$ | $O(H)$ |
| 种群规模调整 | $O(dy\_NP)$ | $O(1)$ |
| **每代总计** | $O(dy\_NP \times D + dy\_NP \log dy\_NP)$ | $O(NP \times D + ANum \times D)$ |
| **总体** | $O(maxFEs \times D)$ | $O(NP \times D)$ |

---

*文档生成时间: 2026-06-17*  
*算法实现: C++ (基于 SHADE 框架独立开发)*  
*测试平台: CEC2017 Benchmark, D=30/50/100*
