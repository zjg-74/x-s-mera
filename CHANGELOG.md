# Changelog

本仓库采用语义化版本。理论条款编号（A1–A12、C1/C2、F1/F2、I1）与《X-S-MERA 修订登记》保持同步。

## v1.6.48 — 2026-10-29

### OPEN-3 实现路线：θ_min/Δ 纳入实例申报表约定声明字段模板（A13 落地）

> 承接 OPEN-3 连续化规格评审定稿（v1.6.31，`32-v3.4-open3` §6）的
> `boundary_map.json` next 字段："实现路线：θ_min/Δ 纳入实例申报表约定声明
> 字段模板（A13 申报规范）；不新增 L1 公理，不改变 A/B 级分级【定理】"。

#### feat: `experiments/a13_first_case/declaration_template.md`（新建）
- 实例申报表**约定声明字段模板**：`open2_declaration`（OPEN-2 A/B）、
  `z_vertex`（R34，约定 1.96）、`eta`（R-042/043/044）、`open3_theta_min`、
  `open3_delta`（OPEN-3，约束 `Δ∈[0,θ_min)`，`Δ=0` 退化二值）。
- OPEN-3 连续化 C 级辨识层三值判读（R34 先例）：`θ_k≥θ_min` 可读 /
  `θ_k∈[θ_min−Δ,θ_min)` 边缘不判 / `θ_k<θ_min−Δ` 闭。
- 诚实边界：Δ 为约定输入、框架不裁定；边缘带不产出闭合结论。

#### feat: `experiments/a13_first_case/a13_contrast.json`（+declaration 块）
- 顶层新增 `declaration` 示范块：`open2_declaration="A"`、`z_vertex=1.96`、
  `eta=0.05`、`open3_theta_min=0.8`、`open3_delta=0.1`，与 OPEN-2/z/η 并列。
  标注"示范性约定声明，非测量断言"。

#### feat: `experiments/a13_first_case/README.md`（+§5 OPEN-3 落地说明）

#### 验证与诚实边界
- `a13_contrast.json` 合法 JSON，`declaration` 块字段完整（open3_theta_min/
  open3_delta/open2_declaration）。
- **诚实声明**：模板为示范性约定声明字段模板，示范块不构成对现实的测量断言；
  连续化不新增 L1 公理、不改变 A/B 级分级【定理】（OPEN-3 规格 §6 定稿）。

## v1.6.47 — 2026-10-29


### B3-ext：LE 内核四类终态机制扩展（L2 仿真层）

> 承接 v1.6.45 B3 Part A 的关键发现：当前 LE 内核（v3.4-le.1）单调稳定，
> 四类终态（稳态/周期振荡/生态锁定/全局灭绝）在内核本身**不可达**（仅扩展
> 原型可及）→ 登记 B3-ext。本轮为内核引入**可选扩展机制**，使四类终态在
> **LE 内核本身**可达。

#### feat: `src/le_coupling_kernel.py`（+ext）
- 新增 `mechanisms` 参数 + `_validate_mechanisms`（`EXT_MECH_KEYS` 校验，非法键
  `ValueError`）。**默认 None 保持 v3.4-le.1 单调稳定内核（向后兼容）**。
- `_life_rhs` 扩展项：
  - **灭绝**：`death_threshold`/`death_rate` → life 低于阈值不可逆加速死亡（吸收态）；
  - **锁定**：`allee_threshold`/`allee_strength` → Allee 双稳（乘 `(1-bio/k)` 有界化）；
  - **振荡**：`oscillation_gain` + `env_r`/`consumer_mortality`/`half_override`/
    `uptake_override`/`conversion` → **纯消费者 R-M 结构**（bio 完全依赖资源摄食
    + 线性死亡，替换 d0）→ Hopf 极限环。
- `_env_rhs` 扩展：振荡机制下资源用 **logistic 再生**（替代线性供给，真 R-M）。

#### feat: `experiments/b3_ext_mechanisms/`（新建）
- `b3_ext_mechanisms.py` 测绘 → `results.json`（seed=20261029，T=2000）。
- **四类代表全部可达**：steady（final 1.309）/ extinction（0.000）/
  lock（0.164）/ oscillation（8.167）。
- **机制分岔扫描**：灭绝阈值 death_threshold≥0.9 → 灭绝；锁定 allee_threshold≥1.0 →
  锁定（低态）；振荡 conversion≥0.8 → 持续极限环（Hopf）。

#### feat: `tests/test_b3_ext_mechanisms.py`（新增 9 项）
- 默认内核向后兼容（单调稳定，末期无波动）；四类可达；三处分岔
  （灭绝/锁定/振荡 Hopf）；非法机制键拒绝。

#### 验证与诚实边界
- `python -m pytest tests/test_b3_ext_mechanisms.py` **9 passed**；全套回归见下。
- **诚实声明**：机制为操作化原型扩展，非对现实的生态学断言；四类标签是数值
  轨迹的操作化分类；**默认内核（无机制）行为不变**。振荡（Hopf）是纯消费者
  R-M 结构的数学性质，不代表认知系统必然振荡。

## v1.6.46 — 2026-10-29


### B-B25 互信息估计器 I(M; Ω | φ)（纯工程闭合）

> GAPS B-B25：跨理解者混合估计偏倚恰等于 `I(M;Ω|φ)`（理论部分已由
> v3.2 审计 S8 闭合——`P := H(M|φ,Ω)` 定义的直接推论），剩余为**纯工程**：
> 互信息估计器。本轮实现估计器并演示。

#### feat: `experiments/b25_mi_estimator/`（新建）
- `mi_estimator.py`：条件互信息 `I(M;Ω|φ) = H(M|φ) − H(M|φ,Ω)` 估计器；
  经验条件熵 + **Miller–Madow 修正** `Ĥ_MM = Ĥ_ML + (K−1)/(2N)`；
  报告分量 `H(M)` / `H(M|φ)`（混合，有偏）/ `H(M|φ,Ω)`（分群体，无偏基准）。
- `b25_mi_estimator.py` 演示（合成标签数据，seed=20261029）→ `results.json`：
  - 场景 A（Ω 影响 M 分布）：`I = +0.3610` → 偏倚可检出；
  - 场景 B（Ω 与 M 无关）：`I = +0.0006` → 可跨群体比较；
  - 场景 C（仅 φ 影响 M、Ω 无关）：`I = −0.0001` → 条件化后偏倚移除。

#### feat: `tests/test_b25_mi_estimator.py`（新增 7 项）
- Miller–Madow 熵（单点=0 / 均匀=log₂K）；条件熵已知分布；
  三场景互信息判定（A>0 / B≈0 / C≈0）；报告分量恒等 `I = H(M|φ)−H(M|φ,Ω)`。

#### 验证与诚实边界
- `python -m pytest tests/test_b25_mi_estimator.py` **7 passed**；
  全套回归 **443 passed**（436+7，无回归）。
- **诚实声明**：估计器用合成标签数据演示与自测，**非真实认知数据**；
  P̂ 的跨群体真实偏倚报告须待 B 线真实样本（外部数据卡点）。
- 推论落地：`P̂` 数值**不可跨群体直接比较**，除非同时报告 `I(M;Ω|φ)`
  （`15-v3.2` §4.3 / `16-v3.2`）。估计器使该偏倚可计算、可报告。

## v1.6.45 — 2026-10-29


### B3 四类终态相图测绘（L2 仿真层）

> GAPS.md B3（归档推论 6：吸引子分类）：体系自发抵达四类终态
> （稳态 / 周期振荡 / 生态锁定 / 全局灭绝）。终态分类是动力学稳定性陈述，
> FOL 谓词层结构性不可达 → 交 L2 仿真层核验。本轮双轨测绘落地。

#### feat: `experiments/b3_attractor_phase/`（新建）
- **Part A · 当前 LE 内核吸引子测绘**：`le_coupling_kernel` v3.4-le.1 耦合网格
  （coupling_gain × d_threshold_c1）36 点扫描 → **全部收敛唯一非零稳态**。
  **关键工程发现（诚实登记）**：LE 内核是**单调稳定系统**（逻辑斯谛增长 +
  Michaelis-Menten 摄食 + 线性资源供给-消耗），**不呈现振荡/灭绝/锁定**——
  B3 的四类终态在当前内核**不可达**，须内核动力学扩展 → 登记为 **B3-ext 待办**。
- **Part B · 扩展生态动力学原型相图**（b3 实验内独立模型，诚实标注"原型"）：
  资源-消费者 LV（logistic 资源自再生）+ 密度约束 + Allee 效应双稳；
  网格（资源再生率 s × 捕食率 a）**四类终态全覆盖**：
  extinction 29 / steady_state 41 / ecological_lock 33 / periodic_oscillation 37。
  相图分区清晰：低 s→灭绝、中等 s→稳态、高 a→振荡（Hopf 极限环）、临界→锁定。
- 产出 `results.json`（`part_A_LE_kernel` + `part_B_eco_prototype` + 分类器默认）。

#### feat: `tests/test_b3_attractor_phase.py`（新增 8 项）
- 分类器四类判别单元测试（steady/extinction/oscillation/lock）。
- 双轨 smoke：Part A LE 单调稳定（全网格稳态）+ Part B 四类全覆盖 + eco 模型有界。

#### 验证与诚实边界
- `python -m pytest tests/test_b3_attractor_phase.py` **8 passed**。
- 诚实声明：Part B 是**操作化原型模型**，非当前 LE 内核；四类标签是对数值轨迹的
  操作化分类，非现实生态学断言；Part A"单调稳定"是对数值网格的观测，未覆盖全参数空间。
- 正式四类终态测绘须在 LE 内核引入振荡/灭绝/锁定机制（**B3-ext**），本实验为其
  可达性提供原型证据。

## v1.6.44 — 2026-10-29

### C4-full 耗散结构完整版形式化（Lean L1）

> GAPS.md 归档推论 4（耗散结构）**完整版**——骨架版（F6→F8、F9）此前可证，
> 完整版"生命以环境熵增为代价维持有序"缺合系统孤立性前提，登记为缺口 C4-full。
> 本轮**采用签名扩展路线**实现（B3 为 L2 仿真层目标，非 Lean，转后续仿真优先级）。

#### feat: `formal/Formal/Basic.lean`
- 新增**二元合系统孤立性**谓词 `IsolatedCouple : Sys → Sys → Prop`（C4-full 前件）。
  **不改单目 `Isolated`**（M1 用），不破坏现有 AttackSurface/推论。

#### feat: `formal/Formal/Axioms.lean`
- 新增**正 Horn 桥接公理** `C4_isolated_entropy_cost`：
  `LifeS → EnvS → IsolatedCouple → (NeedsEnergyForStructure → EntropyNonDec)`。
  无否定前提 ⟹ AttackSurface 第 1 条（正片段无 ⊥ 推导源）保持有效，无需切换模型构造路线。

#### feat: `formal/Formal/Corollaries.lean`
- 更新 C4-full 缺口注释（缺口 → 已实现），新增导出定理 `dissipative_structure_full`。

#### feat: `formal/GAPS.md`
- C4-full 条目更新为**决策（2026-10-29）：采用签名扩展路线**；登记操作化诚实声明。

#### 验证
- `lake build` **通过、无 sorry**（8 jobs，v4.35.0）。
- 操作化诚实声明：Horn 层把"以环境熵增为代价"表达为"合系统孤立 + 生命需持续
  自由能 ⇒ 环境熵不减"——熵"增/不减"的**数值核算**属 L2 仿真层（熵账本），不在 L1 断言数值。

## v1.6.43 — 2026-10-29

### 推论 R1 / R2 落地（宏-槽-事实离散原型）

> 主文 `01-v3.1-main-dichotomy-macro-loading-trichotomy.md` §9 登记的两条工程可验推论
> 首次实现并验证（X2 经核验**不能纯内部推进**，须 P4 加硬=外部数据，未触发）。

#### feat: `experiments/r1_loading_duality/`（NEW）
- **推论 R1**（§4.4 + §9）：加载度 σ(M,F_k) 与工程边界重合度 d 的同构对偶。宏-槽-事实
  离散原型：同一批事实 F 填充相关宏的核心槽 ⟹ 宏间 d=**1.000**、随机配对宏 d=**0.000**
  （均在 F 下计算），separation=**1.000** ⟹ **HOLDS**。
- 诚实边界：离散原型非现实证据；d=1 是"完全填充"结构极限，随机对照证明非恒真。

#### feat: `experiments/r2_failure_separation/`（NEW）
- **推论 R2**（§7.2 + §9）：F1↔发散（¬C1）/ F2↔伪收敛（C1∧¬C2）条件分离。真收敛
  （C1∧C2）、F1 发散（A_k=∅）、F2 伪收敛（收敛到错对象）三态**可分离**，
  conditions_separated=true ⟹ **HOLDS**。
- 诚实边界：离散原型非现实证据；C2 需预设 ground truth。

#### tests
- `tests/test_r1_loading_duality.py`（4 项）+ `tests/test_r2_failure_separation.py`（6 项），
  **10 项**全过（本机 `python -m pytest`）。基线 393 → **403**。

#### docs
- `README.md`：代码对应表新增 R1/R2 两行、运行命令加两条、测试计数 393 → 403。
- `experiments/{r1_loading_duality,r2_failure_separation}/README.md`：各含目的诚实声明、
  模型定义、结果表、诚实边界、运行方式。

#### X2 立项状态（重要更正）
- `formal/X2-cognitive-bridge-charter.md` §7：步骤 2（CognitiveBridge.lean）**未触发**；
  正式立项仍需 **P4 加硬**（多模型、真人 think-aloud/编辑历史、反刍对照）⟹ X2 **依赖
  外部数据**，从"可纯内部推进"清单移出（上一轮误判已更正）。低风险 R5（sign(0)=0）待列。

## v1.6.42 — 2026-10-29

### 认知系统的动力学结构（理论重构 v4.0 单篇）入库

> 用户指定 `workbuddy/2026-10-28-09-09-48/认知系统的动力学结构.md`，将这篇自成体系的
> 独立论文作为 `docs/xs-framework/35-v3.4-cognitive-dynamics-structure.md`（NEW）入库。

#### feat: `docs/xs-framework/35-v3.4-cognitive-dynamics-structure.md`（NEW）
- 《认知系统的动力学结构：一个可测量、可证伪的最小理论》——摘要 + 关键词 +
  **四重门方法纪律**（可表述/可执行/可证明或登记缺口/可判别）+ 七级证据分级。
- 构造闭环（§3–§7）：范畴二分 → 命题 1「约束恰为双」【定理】 → 状态变量
  U := P − Q ≥ 0（P 依赖理解者，跨群体比较需报告 I(M;Ω|φ)）→ 假设 A_D（双速率
  展开-消解，**唯一承重猜想，未证且无已登记证否途径**）→ 规范固定等价于仿射核算律
  dU/dt = ρ_exp·dE/dt + κ_E·ρ_exp·(1−ρ) → 归宿三分/震荡双相/失效 F1/F2 正则分支。
- §8 证据总账表：A_D【猜想】、U≥0【定理】、仿射律条件于 A_D 辖域 U 远离 0、震荡双相
  判别【原型确认】（LLM 18 任务 + S0 真人 3/6，稳健性 75%）、F1/F2 可诊断【定理】
  Lean4 T1–T6 零 sorry、倒 U（R21，Q≈1.0）、快思维本体不可读（R47）、涌现判别【缺口】。
- 作者声明：**未使用六篇散文版系列文本的任何段落**，与散文版共享的仅是最底层范畴直觉，
  全部定义、定律与证据状态以工程修正后的理论重构（v4.0）为准。

#### docs
- `docs/xs-framework/README.md`：文档表新增 35 号行（注明四重门、构造闭环、证据总账、划界）。
- `README.md`：xs-framework 理论文档计数 34 → 35 篇（v3.4 新增 10 篇）。
- 记号/区分核查：本文全部记号（X/S/M/φ/Ω/C_in/P/Q/U/H/I/ρ_exp/ρ_res/κ_E/ρ/s/σ/Ψ/ψ̄/ξ/A_D/
  F1/F2/g/B/E 等）与区分（D-CONV/D-FAIL/D-DRILL/D-SCALE 等）均已登记，**未新增 registry 项**；
  未改任何冻结文档。三道 CI 门本地全绿。

## v1.6.41 — 2026-10-29

### 六篇文章 v2 修订合集入库

> 用户指定 `workbuddy/2026-10-28-09-09-48/六篇文章-v2/`，将六篇理论的 v2 修订版
> 作为 `docs/xs-framework/34-v3.4-six-articles-v2.md`（NEW）入库；对应 20 号
> （八月去重版）的系统升级。凡表述与形式版冲突之处，以形式版为准。

#### feat: `docs/xs-framework/34-v3.4-six-articles-v2.md`（NEW）
- 六篇合集（篇一方法论 V7.0 / 篇二框架本体 V6.0 / 篇三语言 V3.0 / 篇四教学 V2.0 /
  篇五直觉与错觉 V6.0 / 篇六涌现 V6.0），按《六篇文章修订方案》完成 G1–G4 全局修订：
  - **G1-A** X/S 对齐形式版范畴层（L0）二分（X=象对象面 / S=运算运算面）；旧"扩展/固化
    倾向"改写为 L2 层双速率 ρ_exp/ρ_res，不同层不混用。
  - **G2** 熵可测量化：P := H(M|φ,Ω)、Q := I(M;C_in|φ)、U := P−Q ≥ 0，R23 判别器为操作化
    测量程序。
  - **G3** 震荡双相：酝酿型（U 净降，ρ_res>ρ_exp）/ 反刍型（U 净升，ρ_res<ρ_exp）。
  - **G4** 七级证据分级（定理/原型确认/猜想/待验证/缺口/约定/局限）+ 各篇篇首命题地位表。
- 卷首新增"六篇索引与版本"表；正文逐篇含命题地位表、修订说明（V5.0/V6.0→新版）。

#### docs
- `docs/xs-framework/README.md`：文档表新增 34 号行（注明 G1–G4、六篇版本、与形式版对位）。
- `README.md`：xs-framework 理论文档计数 19 → 34 篇（补 v3.4 新增），提及六篇 v2。
- 记号/区分核查：六篇 v2 所用全部记号（P/Q/U/ψ̄/Ψ/σ/ρ_exp/ρ_res/κ_E/A_D/η/s/g/E 等）与
  区分（真/伪收敛 D-CONV、钻研/空洞 D-DRILL、局部/全局 D-LOCAL、外/内失效 D-FAIL、
  同尺度耦合 D-SCALE）均已登记，**未新增 registry 项**；未改任何冻结文档。三道 CI 门本地全绿。

## v1.6.40 — 2026-10-29

### 认知架构内核落地（src/kernel.py）

> 填补仓库登记缺口：README「src/kernel.py ⏳ 计划中」/「认知架构仿真 ❌ 尚不在本仓库」。
> 设计依据：`workbuddy/2026-10-28-09-09-48/认知架构内核-设计草案.md`（v0.1）与
> `kernel_proposal/kernel.py`（提案稿）。实现设计草案 §6「第一步只做这些」第 1–4 条。

#### feat: `src/kernel.py`（NEW）
- **内核骨架**：状态 `(E, U)` + 规范固定参数 `GaugeFixedParams(s, rho_exp, rho, kappa_E)`
  + 仿射积分器 `step`（v3.2 主定理 dU = rho_exp·dE + κ_E·rho_exp·(1−ρ)·dt）+ A_D⁺ 的
  `g(U)` 结构约束（g_saturated/g_linear，B-B43 处置：不实现形式拟合）。
- **内核级自证伪**：`t_affine`（纯 A_D 共线性残差 ≤2e-15，R30 基准）、
  `t_gauge`（G_λ 规范不变性，定理 6.12）。
- **D-SCALE 判别接口迁入**：`fit_beta_with_se` + `dscale_verdict`（三值：
  mid_drill / mid_hollow / vertex_band，R34 禁区 {|β̂|≤z·SE}，z=1.96 约定性）+ `dscale_classify`。
- **R23 输出适配器**：`to_r23_json`（内核轨迹→analyze_trajectories.py 格式；
  权重透传、不虚构）+ `rho_res_verdict`（复刻 R23 的 ρ_res>ρ_exp 走向分类）。
- 诚实边界（写死进模块 docstring）：数值原型非经验验证；A_D 未证（B-B22）+ 单一驱动
  附款（S2，调用方契约）；快思维本体不可测（R47/Q0）；g(U) 形式是约定（R45）；读数
  条件于 A_D 与声明的尺度。

#### tests
- `tests/test_kernel.py`（NEW，14 项）：参数校验、step 推进与记录、T-affine 纯 A_D 共线
  残差≈0、T-gauge 一致/偏差、D-SCALE 三值判读 + 禁区随 z 变宽、楔形方向（非共线观测
  E=k²）、R23 适配器格式与权重透传、ρ_res 走向分类与无震荡不可判。仅标准库、调底层
  函数、零文件写入。

#### docs
- `README.md`：顶部版本 v1.6.39 → **v1.6.40**；Repository state 把「认知架构仿真 ❌ 尚不在
  本仓库」更新为已落地，增补 kernel.py 行（含诚实边界）。
- 未改任何登记件/冻结文档（内核所用记号 E/U/ρ_exp/ρ/κ_E/s/g 均已登记）；三道 CI 门本地全绿。

## v1.6.39 — 2026-10-29

### 补录 R46 / R47 运行记录（证据补录，无代码/结论变更）

> 仓库更新检查：工作树仅 3 个未跟踪文件。其中 `r23_real_cognitive/r8b_uniqueness.json`
> 与已入库的 `r8_gain_stepchain/r8b_uniqueness.json` **内容完全相同**（重复副本，不入库）；
> 其余两个为 R46 / R47 的运行 stdout 记录，补录入库作为运行证据。
> 远程 main 与本地 HEAD 一致（= a7b7d3cb，v1.6.38），无远程新提交。

- `experiments/r23_real_cognitive/r46_run.txt`（NEW）：R46 语言化序列化算子运行输出。
- `experiments/r23_real_cognitive/r47_run.txt`（NEW）：R47 δ-时间化检验运行输出。
- R46（v1.6.35）、R47（v1.6.36）的脚本 / report / json 均已在此前入库，本段仅补录运行记录。

## v1.6.38 — 2026-10-28

### 纯内部批次：边界图机检接入（A10 残留）+ B-B32 震荡判据一致性检验

> 承 v1.6.37 用户"继续清纯内部项"。两项均纯内部可解（S2 / 机检），真实数据仍待 B 线采集。

#### feat: A10 残留落地 · 边界地图机检（`src/xs_boundary_map_check.py` NEW）
- 让边界地图 `docs/xs-framework/boundary_map.json` 从"仅登记"推进到"被机器看住"
  （对应 27 §6.5 A10 残留、逻辑桥"让地图被机器看住"）。
- 校验器（纯标准库，仿 xs_formal_checks 模式）：顶层必填字段、dimensions/r041_table/
  open_cells/gap_index 必填与 id 唯一、gap_index 已解决格 resolved_in 指向仓库内真实文件、
  状态词软校验；error>0 退出码 1。
- 本地运行 0 错误 0 警告（当前登记件完整）；`boundary_map.json` 升级 **v0.6**（purpose 更新为
  "机检已接入"，revised 追加 2026-10-28）。
- 接入 CI 门 `xs-gates.yml` 第四步（`python src/xs_boundary_map_check.py`）。
- 测试：`tests/test_xs_boundary_map_check.py`（11 项，纯标准库、调底层函数、零文件写入）。

#### feat: B-B32 震荡两种判据一致性检验（`experiments/b32_shape_prob_consistency/` NEW）
- 缺口：补篇七 B-B32——推论 6.11 形状判据（震荡 = 局部收敛被重开，重开次数可数）vs
  主文概率判据（U>0 且仍在下降），一致性未检验。
- S2 合成候选竞争模型（固定种子 20261028，纯标准库）：网格 20 条轨迹（c×o），
  对每条分别用两判据判"震荡/非震荡"，比较一致性与不一致形态。
- **结果：一致率 0.75，双向均有反例，一致性不成立**——
  - 反复「局部收敛被重开」的经典震荡中，U 被重开反复抬高、不满足"U 仍在下降"（概率判据漏判）；
  - 从未到 m=1 的未决振荡中，无局部收敛可被重开（形状判据漏判）。
- 结论：形状判据锁定"局部收敛被推翻"的结构事件（要求先到 m=1），概率判据锁定
  "未决熵收敛趋势"（要求 U 下降）。**二者互补而非冗余，需合并覆盖**；主文若以"震荡"为统一
  构念需声明主判据。
- 输出 `b32_consistency.json` + `report.md`；08 supp7 缺口表 B-B32 更新为"已检验"。

#### docs
- `CHANGELOG.md`：本段（v1.6.38）。
- `docs/xs-framework/boundary_map.json`：v0.6（机检接入登记）。
- `docs/xs-framework/08-v3.1-supp7-micro-reversal-locality.md`：B-B32 状态 → 已检验。
- `README.md`：顶部 v1.6.37 → **v1.6.38**，Repository state 增补两行。

## v1.6.37 — 2026-10-28

### B-B9 唯一性推进：势消除混沌的唯一机制归约（M1/M2/M3 消融判别）

> B-B9（补篇一缺口表）：未证明 M1 是势消除混沌的**唯一**机制，M2/M3 未被排除。
> R8（v1.6.27 起，`r8_gain_stepchain`）只证了 M1 的**存在性**（周期倍增级联 V1–V4）。
> 本版在 R8 同一原型上做**机制消融判别**，把"唯一性"推到位移。

#### feat: `experiments/r8_gain_stepchain/r8b_uniqueness.py`（NEW）
- 宿主：R8 同一 logistic 型离散原型 + 对称双候选扩展（g1=g(1+s), g2=g(1−s)）。
- 四判据（纯标准库、固定种子 20261028）：U1 M1 充分且合规 / U2 非法通道排除 /
  U3 M2 合法形式归约 M1 / U4 M3 正交（只压缩寿命）。
- 输出 `r8b_uniqueness.json`；`report.md` 追加「R8-B9」小节。

#### 结果（四判据全真 ⟹ 框架内唯一性成立）
- **U1**：M1 充分合规——混沌窗内增益调制把 λ 由正转负（g=3.8 λ=0.43 → g=3.5 λ=−0.87），
  g 仅乘性（A11-R 合规）。
- **U2**：非法通道排除——势放进状态导数（加性偏向 b）→ A11-R 合规标记 FAIL，通道被排除。
- **U3**：M2 合法对称破缺 = 增益不对称 = M1 推广——s 只进乘性增益、单调驱动动力学
  （λ(s) 0.42→1.51，无加性/梯度项）⟹ M2 ⊆ M1。
- **U4**：M3 正交——强制归宿不移动 g_chaos（恒 3.57），只单调压缩混沌寿命
  （300000→51）。

#### 结论（诚实）
在 **A11-R（势-流分离）约束下**，唯一能改变"混沌是否出现（拓扑）"的势机制是 **M1**；
M2 合法形式归约 M1、非法形式被排除；M3 与拓扑正交（管"持续多久"而非"是否出现"）。
**B-B9 由"完全开放"推到位移**：框架内唯一性已证；剩余"跨框架全局唯一性"非本框架可裁，
转登记为外部开放项（不越界声称全局唯一）。

#### docs
- `docs/xs-framework/27-v3.4-boundary-map.md`：B-B9 登记由"仍开放"更新为"框架内唯一性已证
  （R8-B9 v1.6.37），跨框架全局唯一性转外部"（§5 开放格 / §6 逻辑桥 / 下一颗钉子）。
- `experiments/r8_gain_stepchain/report.md`：追加 R8-B9 小节。
- `README.md`：顶部 v1.6.36 → **v1.6.37**，Repository state 增补 R8-B9 行。

## v1.6.36 — 2026-10-28

### R47：δ-时间化检验——震荡/方向是否时间化算子的产物（B 线 P4 方法论）

> 承接用户命题（2026-10-28）：「快思维 A→B 是瞬时跃迁（τ=0、δ、无时间结构）；震荡不是
> 认知本体，而是时间出现在函数中才产生的表观；时间化了就不是快思维，而是慢思维」。
> R46 检验了"语言化是否翻转/抹平方向"；R47 更根本——从**无内部时间结构的 δ 事件**出发，
> 检验"震荡/方向是否完全是时间化算子的产物"。

#### feat: `experiments/r23_real_cognitive/r47_delta_temporalization.py`（NEW）
- 快思维 δ 事件：瞬时选择（τ=0、无内部步），只携带候选集中度结构（主导占比 p 演化）。
  三类：酝酿（p: 0.25→0.85，熵降）、反刍（0.85→0.25，熵升）、随机（恒定，无结构）。
- 时间化算子：δ 事件序列展开成时间轨迹（每事件 m 子步、焦点沿主导类演化、focus 控制
  锐化）；总质量 T=100（避免量化平台/熵阶梯跳变）。
- 三问检验 Q0/Q1/Q2 + `r47_delta_temporalization.json` / `r47_report.md`。

#### 结果（三层，命题被精确裁决）
- **Q0 必要条件**：纯 δ 单点（无时间化）判别 **30/30 无震荡段、无判定** ⟹ 判别器需时间
  结构，**快思维本体（瞬时选择）不可测**——实证"时间化了就不是快思维，而是慢思维"。
- **Q1 不造方向**：随机 δ + 时间化方向检出率**全 0**（判空转为主）⟹ 时间化**不凭空造方向**，
  否定"震荡/方向纯是时间化伪影"的最强形式。
- **Q2 结构恢复**：酝酿/反刍 δ 跨算子（m=1–5, focus≤0.6）方向保持率**全 1.0**；仅极端全聚焦
  （focus=0.9，单候选 r=1）崩到 0。

#### 结论（诚实）
判别器测的方向**部分反映慢思维真实结构、非纯时间化伪影**；但快思维本体方向不可测，
读数仍登记为"**慢思维（时间化认知）轨迹**"。采集约束与 R46 一致：引导"保持多候选并存地
出声"、规避单线聚焦；触达快思维需非时间化通道（反应时/生理/鼠标速度/隐式启动）。

#### docs
- `experiments/r23_real_cognitive/README.md`：追加「R47 δ-时间化检验」小节。
- `README.md`：顶部 v1.6.35 → **v1.6.36**。

## v1.6.35 — 2026-10-28

### R46：语言化序列化算子 + 判别方向稳定性（B 线 P4 方法论，承接用户批评）

> 用户批评（2026-10-27）：「采集录音得到的也是污染了的数据：认知一旦被语言化，认知过程
> 就是一个慢思维的显化」。本版把该批评转成可证伪的 S2 检验：在已知方向（酝酿 U 降／反刍
> U 升）的合成轨迹上施加**语言化算子**（时间细分 m + 焦点放大 focus），对比判别器
> （与 `analyze_trajectories.py` 同逻辑）读出的方向。

#### feat: `experiments/r23_real_cognitive/r46_language_serialization.py`（NEW）
- 合成已知方向轨迹：4 候选类、8 步、伪计数总量 10；酝酿（主导权重 0.25→0.85）与反刍
  （0.85→0.25）各 12 条，原始判别全对。
- 语言化算子两模式分离：**模式A·纯时间细分** `subdivide_keep`（等比例拆 m 子步、不改分布，
  对照组）＋**模式B·聚焦锯齿** `subdivide_focus`（子步向焦点类放大 boost=1+focus、保留原始
  结构）。参数扫描 m∈{1,2,3,5}、focus∈{0,0.3,0.6,0.9}。
- `r46_language_serialization.json`（meta / original_verdict / lang_verdict / verdict /
  conclusion）/ `r46_report.md`。

#### 结果（✅ 方向稳健，污染=抹平而非翻转）
- **翻转率全部 = 0**：语言化在任何参数下都不把酝酿判成反刍（或反之）。
- **现实采集参数（m=2/3, focus=0~0.3）方向保持 100%**（24/24）。
- **抹平（flat/空转）仅发生在病理性语言化**：高聚焦（focus≥0.6 且 m≥3）或碎片化（m≥5，
  每子步伪计数太少→无震荡段）。

#### 机制含义（对 B 线 P4 / think-aloud 方法论）
1. think-aloud 判**方向**（酝酿/反刍）对语言化显化稳健——污染不改变方向，只可能抹掉方向。
2. 污染的真实形式是**抹平→空转**：被试**单线聚焦**（只讲当前一条路、忽略仍在考虑的多条）
   或语言**碎片化**时，判别退化为无震荡段/空转。
3. **采集协议可操作约束**：引导语应**强调「保持多候选并存地出声」**——否定 A、尝试 B 时，
   仍说出 C 等仍在考虑的路径，避免单线聚焦（与 `thinkaloud-math-script.md` 的"多路径/多设元"
   设计一致，现获量化依据）。
4. 方向稳健 ≠ 读数精确：幅度/τ\* 仍失真（语言化改变时间尺度），需非语言化通道或独立标定。

#### docs
- `experiments/r23_real_cognitive/README.md`：追加「R46 语言化序列化算子」小节。
- `README.md`：顶部 v1.6.34/393 → **v1.6.35**。

## v1.6.34 — 2026-10-27

### R8 / R11 / R44 结构测试补齐 + 台账同步 + 两处既有校验缺陷修复

> 承接 v1.6.33。三项早于 v1.6.31 落地的实验（R8/R11/R44）此前无 pytest 覆盖，
> 且台账停在 v1.6.21。本版补齐测试、同步台账，并**修复两处既有校验失败**（非本版测试引入）。

#### test: 新增三组结构测试（共 15 项，仅标准库，锁结构不锁数字）
- `tests/test_r8_period_doubling.py`（**8 项**）：锁 logistic 形式（g 仅乘性）、迭代值有界
  ∈[0,1)、g=4.0 必混沌（λ>0 且无非平凡短周期）、g=2.5 周期 1、`_coarse_scan()` 复现 main
  扫描逻辑（不调 main，避免写文件污染仓库）、V1 周期倍增链 1→2→4→8、V2 混沌带宽≥20、
  V3 g_chaos∈[3.4,3.7]、M1 机制 V1∧V2∧V3∧V4。
- `tests/test_r11_pseudo_drift.py`（**5 项**）：调 `run_arm_g(arm,g)`（不调 main）。锁两臂均
  有 F2 字段、P1 legal F2 极差≤4pp、P2 violate F2 峰值−legal 峰值≥20pp 且自基线升起、
  violate 每档 g 的 F2 均高于 legal 峰值且上升幅度 >2×legal 波动 ⟹ 推论 4.1 成立（越权
  ⟹ F2/ghost 伪力）。
- `tests/test_r44_eta_principle.py`（**2 项**）：调底层 `make_trajectory/ratio_curve/tau_star`。
  锁 τ* 随 η 单调不增、三场景在落带区间内存在使落带率≥0.90 的 η（可靠子区间非空）。

#### fix: 两处既有校验失败（ regression 回归时发现，均由更早版本引入）
1. **S006（记号表）**——`symbol_registry.json` 中 `eta` 的 `derives_from` 在 v1.6.32 被填入
   4 条**描述性字符串**（可接受错误率语义／场景落带区间／可靠子区间／判定灵敏度），效果是
   ① 引用未登记记号触发 S006；② 违反「η 是约定性输入、`derives_from` 必须为空」契约
   （`test_a7_core_min.py::test_eta_is_registered_as_convention_input` 对此有显式断言）。
   **处置**：`derives_from` 清空为 `[]`，四条外部标定来源（R-042 带内不敏感／R-043 落带区间／
   R44 可靠子区间／OPEN-2 观测协议约定）**移入 `definition` 字段显式登记**（U3：C 级量须登记
   其外部标定来源），信息未丢失。
2. **G002（冻结窗）**——`docs/xs-framework/06-v3.1-supp5-bridge-constraint-strength.md` 在
   **v1.6.28（`089cf90`）** 被修订（§自审第 3 条：R21 复核 σ 离散化量子 1/8，N_FEAT=16／
   N_MACRO=40 下倒 U 保持、峰值仍在 Q≈1.0，已排除伪像），但**未重新冻结**，manifest 仍锁
   v1.6.12 哈希 ⟹ 冻结窗被破坏。已核查 git diff 干净（工作区未被改动）、该修订为有意且
   唯一（4faf5ea→089cf90 之间仅此一处）⟹ **补重冻结**：哈希刷新至 v1.6.28 内容
   `e13b59d1481a…`，并在 `reason` 登记 v1.6.28 修订说明。

#### docs: 台账同步
- `CITATION.cff`：1.6.21 → **1.6.33**（date-released 2026-10-27）。
- `README.md`：顶部 v1.6.21/368 → **v1.6.34/393**。
- `boundary_map.json` → **v0.5**：修正 `version` 字段仍为 "0.4" 的自相矛盾（文件内已写 v0.5
  内容但 version 未同步）；OPEN-4 状态「待验证」→「已闭合（v1.6.22）」；sufficiency 改写
  （CI 门已在 GitHub 启用，不再写"尚未启用"）。
- `27-v3.4-boundary-map.md` → **v0.5**：标题/版本声明 + §3.4③、§4.2、§4.5、§4.6 + §5 五处表格行
  + §6 两条 + §5 总评段整段重写（R-043 已解决／R-044 锚定／OPEN-3 定稿／B-B44 闭合／
  B-B43 降级／B6 三类探索已启动；剩余开放格收敛为 B-B9 唯一性、X2 触发、CI 持续验证）。
- `experiments/r23_real_cognitive/thinkaloud-math-script.md`（NEW，并入）：B 线 P4 M 系列
  数学 think-aloud 采集脚本（判别信息三要素 + 10 题 + 归约规则 + 8/9 触发线）。

#### 回归
- 全套（除 R29）**378 passed**；R29 单跑 **15 passed** ⟹ **393 passed**（含 numpy 环境口径）。
- 两道形式校验器 `xs_formal_checks.py` / `xs_trace_gate.py` **均 exit=0**。

## v1.6.33 — 2026-10-27

### B-B43（g(U) 形式）降级：R45 真实数据可行性评估（加噪稳健性 + 形式辨识）

> 承接 v1.6.32。用户"B-B43"：评估 R32 反演能否由真实 (E,U) 数据闭合 g(U) 的具体形式。
> 纯内部可解（合成数据扩展 b3_residual_inversion）。

#### feat: `experiments/b3_residual_inversion/r45_noise_formid.py`（NEW）
- 复用 b3 内核，加 **U 观测噪声**（σ∈{0,0.01,0.02,0.05,0.1}，M=40）+ **去噪平滑分支**
  （残差滑动平均 w=11）+ **形式辨识**（反演形状 vs 候选族 {linear,square,tanh,saturated}）。
- `r45_noise_formid.json` / `r45_report.md`。

#### 结果（诚实负面）
- **加噪稳健性**：裸反演 σ=1% 即崩（corr 1.0→0.05，差分放大噪声 50×）；滑动平均无法挽救
  （σ=0.01 仅 0.26–0.37）；容忍 σ 上限 = 0.0。
- **形式辨识**：无噪识别正确但 margin≈0（形式族在 U 窄范围形状高度相关）；有噪识别率 50%。

#### 结论
- **B-B43 降级**：g(U) 具体函数形式**不能由现实 (E,U) 数据可靠辨识** ⟹ 正解 = 理论约束
  （结构 g(0)=0、单调、渐近 g(U)→1 + 变异形状 R32/R43），非实证形式拟合。修正 R32
  "套用真实数据即可"的乐观预期。
- **D-SCALE 不受影响**：g 不进判定流程（18 章 §4.4.1）。
- **辖域限定是本质的**：补篇九"U 远离 0"限定不能靠实测 g 放松，必须靠结构约束维持。

#### docs
- `symbol_registry.json`（更新）：g 定义补渐近约束 + R45 不可辨识结论。
- `boundary_map.json`（更新）：B3 status 补 B-B43 降级。
- `README.md`（更新）：缺口表 B-B43 降级为"结构约束+变异形状"。
- `CHANGELOG.md`（本条目）。

## v1.6.32 — 2026-10-27

### B-B44（η 取值原则）锚定闭合：R44 可靠子区间 + 取值原则定稿

> 承接 v1.6.31。用户"B-B44"：给出 η（尺度分离阈值）的取值原则，把"取值原则开放"
> 推进为"已锚定"。纯内部可解（复用 R-042/R-043 内核 + 噪声注入 + Bootstrap）。

#### feat: `experiments/r44_eta_principle/`（NEW）
- `r44_eta_principle.py`（纯标准库，固定种子 20261027）：复用 R-042/R-043 内核，加噪声
  （σ∈{0.05,0.10}）+ Bootstrap（M=200 轨迹/场景/σ）估计落带率 P(τ*∈(8,128))，在落带区间
  （R-043）内压出**落带率≥90% 的可靠子区间**。
- `r44_eta_principle.json` / `report.md`（诚实边界：可靠子区间为原型内结果，REL_TH=0.90
  为约定水平；D-SCALE 仍 provisional，B-B43 未解）。

#### 结果
| 场景 | 可靠子区间（落带率≥90%） |
|---|---|
| A 同尺度 | [0.031, 0.046]（≈落带区间，最差落带率 0.975/0.920） |
| B 跨尺度弱 | [0.032, 0.071]（**下界收窄**；近 0.031 处落带率仅 0.61–0.72） |
| C 结构性 | [0.061, 0.086]（=落带区间，最差落带率 1.000） |

#### feat: `docs/xs-framework/33-v3.4-eta-scaling-principle.md`（NEW）
- **η 取值原则 P1–P5 定稿**：① 语义（可接受错误率，不投射外部显著性水平，守 forbidden）；
  ② 落带（R-043）；③ 可靠子区间（R44，落带率≥90%）；④ 带内不敏感（R-042）；⑤ 多场景
  （交集为空 → 场景相关声明，全局单值取保守 A∩B + 显式声明牺牲场景）。

#### docs
- `symbol_registry.json`（更新）：η definition 补取值原则，derives_from 补四依赖。
- `trace_gaps.json`（更新）：η severity "中"→"已锚定（约定根）"，reason/todo 更新
  （η 仍为约定根，取值原则闭合，不删除）。
- `reef_shturl.json`（更新）：R-043 status 追加 R44；新增 **R-044** 条目。
- `boundary_map.json`（更新）：B3 尺度边界 status 补 R-043 解决 + R44 锚定。
- `README.md`（更新）：v1.6.32 事件登记；文档表增 `33`；代码表增 R44 行。
- `CHANGELOG.md`（本条目）。

## v1.6.31 — 2026-10-27

### 纯内部批次：OPEN-3 规格评审 + R8/R11 判决性实验（M1 机制 / 推论 4.1）

> 承接 v1.6.28 纯内部可推进项。用户"先推其他纯内部项"：OPEN-3 三待评审点拍板定稿，
> R8（增益-步长）与 R11（增益越权）两条判决性实验在离散迭代原型上直接跑通。

#### feat: `docs/xs-framework/`
- `32-v3.4-open3-continuous-observability-spec.md`（更新）：追加 **§6 评审结论 v1.6.31**
  （点1 θ_min 并入实例申报约定字段；点2 弱可观测 → C 级三值降级，R34 先例；点3 程度化
  严格限 C 级辨识层，不连带 B 级 z 容差带）。OPEN-3 由"待评审"→"已评审定稿"。
- `boundary_map.json`（更新）：revised **v0.5**、OPEN-3 状态、sufficiency verdict
  （OPEN-4 未验证 + CI 门未启用仍拦"解锁"）。

#### feat: `experiments/r8_gain_stepchain/`（NEW）
- `r8_period_doubling.py`（纯标准库）：logistic 型最小离散原型，语义势增益 g 扫过
  [2.5,4.0]×1001，周期检测 + Lyapunov。**完整周期倍增级联 1→2→4→8→16→32→混沌**；
  分岔点 2.9965/3.448/3.5455/3.565/3.5695 vs logistic 理论 3.0/3.449/3.544/3.564/3.5699，
  混沌阈值 3.571（理论 3.5699）。V1–V4 四判据全真 ⟹ **M1 机制工程验证**。
- `r8_results.json` / `r8_bifurcation_tail.json` / `report.md`（诚实边界：最小原型≠r18
  完整宏加载原型；阈值网格离散化误差 <0.03%；不声称认知必然混沌）。

#### feat: `experiments/r11_pseudo_drift/`（NEW）
- `r11_pseudo_drift.py`（纯标准库）：复用 r18 宏加载原型，双臂对照——**合法注入**
  （g 只调幅度、不进状态导数）vs **越权注入**（g 进状态导数，冒充漂移）。g∈{1..5}、
  800 试次/格、固定种子 20261027。**越权 F2 峰值 66.1%（g=3）vs 合法 12.9%（+53.2pp）**；
  合法近似平稳（极差 2.1pp）。P1/P2 全真 ⟹ **推论 4.1 成立（越权 → F2/ghost 伪力）**。
- `r11_results.json` / `report.md`（诚实登记：g=4/5 回落 = 冒充漂移饱和次级效应）。

#### docs
- `README.md`（更新）：v1.6.31 事件登记；代码层对应表增 R8/R11 两行；R18 建议段改"已执行"；
  边界地图 OPEN-3 状态推进。**B-B9 部分闭合**（M1 存在性证，唯一性仍开放）。
- `CHANGELOG.md`（本条目）。

## v1.6.30 — 2026-10-27

### P4 B 线加硬 · 真人 think-aloud 数学 M 系列样本（batch1/2，7 题）

> 承接 v1.6.29 单题后，用户按《数学 think-aloud 采集脚本》录 6 段（M1/M2/M3/M4/M5/M7，
> 10月27日 3–8.m4a）。归约 7 段 M 系列 think-aloud 并判别。

#### feat: `experiments/r23_real_cognitive/`
- `raw_transcript_human_1027_m{1,2,3,4,5,7}.json`（6 个，NEW）：M1–M7 真人出声思维 ASR 转写。
- `trajectory_human_1027_m1m2.json` / `trajectory_human_1027_m3457.json`（NEW）：归约轨迹
  （候选=解题路径，每步 note）。
- `results_human_1027_m1m2.json` / `results_human_1027_m3457.json`（NEW）：R23 判别输出。
- `README.md`（更新）：新增"S0 think-aloud 数学 M 系列样本"章节 + 诚实边界。

#### 判别
- 可判 **7：反刍 4 / 酝酿 2 / 空转 1**。
- **关键发现**：方向几乎完全由"被试是否解出"决定——未解出（M1/M2/M4）→ 反刍（发散）；
  解出/收敛（M3/M5/M7）→ 酝酿/空转。验证"未解出偏差"：数学题上"解题成败"成为主混淆变量。
- 对 X2：纯方向反刍略多但受未解出主导，非干净信号；需严格单人 think-aloud + 控制难度
  分离"酝酿/反刍"本质与"解题成败"。
- 诚实边界：5/7/8 多说话人（师生/同伴对话式）非严格单人；M5 结论错误（81 vs 9）仍判酝酿
  ——判别读 U 走向、非答案正确性。

## v1.6.29 — 2026-10-27

### P4 B 线加硬 · 真人 think-aloud 单题样本（10月27日.m4a）

> 承接 v1.6.26 真人6题（10月26日.m4a）后，用户提供**首段严格 think-aloud** 录音
> （10月27日.m4a：数学题出声思维，151.7s，单说话人）。归约为 S0 think-aloud 样本并判别。

#### feat: `experiments/r23_real_cognitive/`（1027 think-aloud）
- `raw_transcript_human_1027.json`（NEW）：真人 think-aloud ASR 完整转写（含 speaker 时间戳，mediakit）。
- `trajectory_human_1027.json`（NEW）：单题归约轨迹（候选=解题路径，每步 note）。
- `results_human_1027.json`（NEW）：R23 判别输出。
- `README.md`（更新）：新增"S0 think-aloud 单题样本"章节 + 诚实边界。

#### 判别
- 题目：漏写两个三位数间乘号得六位数，恰为两数乘积的三倍。**首段严格 think-aloud**（非日常
  快速作答），含真实候选震荡（设六未知数 abc/def →"太未知数太多"→改设 a、b；策略
  方程→试凑→整除多路并行）。
- 结果：最长震荡段 k3–k8（6 步），U 净 **+0.6144**（升）→ **ρ_res<ρ_exp（反刍）**。
- 诚实边界：单题单样本不能单独定论；归约粒度（候选=解题路径）影响结果；开头设元震荡
  （k0–k2，U 降）段较短未被取为最长段。价值：判别信息高于课堂板演，作反刍侧对照，
  需更多 think-aloud 样本形成多数方向。

## v1.6.28 — 2026-10-27

### 纯内部可推进项批次（OPEN-3 规格 / σ 复核 / 事件公理 / R23 稳健性 / 25-audit）

> 承接 GAPS 形式化（v1.6.27）后，按次序推进全部"不需外部数据"的纯内部项。共 5 项，
> 全部完成或登记实际状态。

#### feat: `docs/xs-framework/32-v3.4-open3-continuous-observability-spec.md`（NEW，OPEN-3 规格）
- OPEN-3（B1 维"开/闭二分理想化"）连续化规格：以可观测度 `θ_k ∈ [0,1]` 替代二值 `b_k ∈ {0,1}`；
  可读判据连续化为 `∃k, θ_k ≥ θ_min`，θ_min 为观测协议声明阈值（与 z/η/OPEN-2 同类）。
- 论证：连续化**不改变 A/B 级分级**【定理】（A/B 判读不依赖 s(k) 程度）；只给 C 级辨识层加
  约定阈值参数，不新增 L1 公理。留 3 个待评审点（θ_min 并入实例申报字段 / 弱可观测降级规则 /
  程度化作用范围）。
- `boundary_map.json`：OPEN-3 状态 → "规格已拟（v1.6.28，32-v3.4-open3），待评审"。

#### σ 量子化复核（06-supp5 §7 待办回写登记已解决）
- 确认该待办**已被 R21 §3.1（`r21_highres.py`，N_FEAT=16/N_MACRO=40）覆盖**：倒 U 保持、峰值
  仍 Q≈1.0，排除 σ 离散化量子伪像。06-supp5 §7"未排查"标记回写为"已复核"，总账同步（复用
  既有数据，不重做）。

#### 事件公理定义（23-audit §九）
- 选型 (b) 有限 Boole 代数：「事件公理」= 事件空间原子分解基底（n 个互斥穷尽原子 → 2ⁿ 事件），
  n 由 Stone 表示唯一性内禀；组合并/交/补、系数仅 0/1、不用加权标量。依据：事件自带 0/1 语义
  排除 (a)；n 内禀排除 (c)。

#### R23 判别器稳健性（`r23_robustness.py` + `r23_robustness.json`）
- Bootstrap 权重扰动（固定种子，B=200）：全局多数方向"ρ_res>ρ_exp" **75%** 保持（未达 >90%
  强一致，属中等稳健）；部分任务（T7/R6/R9）扰动下翻转。诚实登记：判别器方向基本稳定、
  非全部任务鲁棒。

#### 25-audit 待办登记实际状态
- #3 补引用、#4 知识注入测试、#5 缺口 1：确认**均已由后续文档完成**（`12-meta-rules` 文献收录、
  `10-v3.0` 主文档 P2 表 + 判定段 + 缺口 1 登记）。#1 反馈协议（Hesse 语气）已内嵌 §四。
  #2（观察者确认辖域）仍需外部观察者，留待。

## v1.6.27 — 2026-10-27

### GAPS B1/B2 桥接公理论证 + AttackSurface 一致性检查（纯内部形式化）

> 用户指定推进纯内部形式化缺口闭合。对 GAPS.md 的 B1（衰老）/B2（演化）缺口做桥接
> 公理决策并编码，重跑 AttackSurface 一致性检查。

#### 决策（论证结论）
- **B1（衰老）→ 降级**，不引入桥接公理（B1a/B1b 均不采纳）。
  理由：必然老化是**无资格条件的全称动力学必然**，强制入 L1 引入未经验证的强断言，
  应留 L2 仿真层（与 B3/X1 同构）。降级形态 = `aging_mechanism_exists`（F4 重述，已编码）。
- **B2（演化）→ 条件化引入**。与 B1 不同，定向演化有明确资格条件（R5：Ns≳5），
  可显式化为正前提。新增 `DirectedEvolution` + `EvoQualified`（C_q 资格抽象）谓词，
  正 Horn 桥接公理 `B2_evolution`，资格显式化与 R5 一致。
- **对称性洞察**：B1 降级 vs B2 引入的差别在于**资格条件可否显式化**——B2 可，B1 不可
  （无条件全称必然）。

#### feat: `formal/`
- `Basic.lean`：新增谓词 `DirectedEvolution`、`EvoQualified`（新概念申请，见 GAPS B2）。
- `Axioms.lean`：新增正 Horn 桥接公理 `B2_evolution`。
- `Corollaries.lean`：登记 B1 降级；新增定理 `directed_evolution`（由 B2_evolution 导出）。
- `AttackSurface.lean`：登记重跑记录——B2 为正 Horn（无否定前提），第 1 条一致性保持有效。

#### 验证
- `lake build` 通过（Lean v4.35.0），全工程无 sorry/admit（G1-a 判据满足）。
- AttackSurface 攻击面分析：正片段无 ⊥ 推导源，无需切换到模型构造路线。

## v1.6.26 — 2026-10-27

### P4 B 线加硬 · 真人录音 S0 样本（单人6题判断）入库

> 承接 v1.6.25 的 P4 外部加硬主线程。用户提供真人录音 `10月26日.m4a`（单人针对
> 6 道生活化/思辨题逐一输出判断与选择逻辑，172.5s，3 说话人）。

#### feat: `experiments/r23_real_cognitive/`（真人 S0 数据）

- `raw_transcript_human_6q.json`（NEW）：真人录音 mediakit ASR 完整转写（含 speaker 时间戳）。
- `trajectory_human_6q.json`（NEW）：6 题诚实归约轨迹（择校=真实 A/B 权衡；消费=条件分叉；
  时间分配=排序；悖论/沙堆/平行线=单步直接作答）。
- `results_human_6q.json`（NEW）：R23 判别输出。
- `README.md`（更新）：新增"S0 真人录音样本"章节。

#### 判别结果
| 题 | 震荡段 | U 净变 | 判定 |
|---|---|---|---|
| H6q-1 择校 | [0..3] | −0.425 | ρ_res>ρ_exp |
| H6q-2 消费 | [0..2] | −0.049 | ρ_res>ρ_exp |
| H6q-3 时间分配 | [0..2] | +0.055 | ρ_res<ρ_exp |
| H6q-4/5/6 | — | — | 无震荡段 |

可判 3/6：酝酿 2 / 反刍 1。

#### 诚实边界（关键）
- **仅 3/6 题构成可判震荡段**；题4/5/6 单步直接作答（"有可能真话"/"主观判断"/"不会做"），
  无多候选并存，判别返回无判定——真实反映日常作答多为**快速单线决策**。
- 唯一**真实 A/B 权衡**是择校题（判酝酿）；消费/时间分配为条件分叉/排序，权重主观、判定弱。
- 采集为**非严格 think-aloud**。
- ⟹ S0 真人数据补充但**判别力有限，不足以独立触发 P4 外部加硬/正式立项**；
  与课堂板演结论一致——可靠 S0 判别仍需严格 think-aloud 或专门多解发散题。

## v1.6.25 — 2026-10-25

### P4 B 线加硬 · 课堂视频 S0 素材 + 采集协议 + 板演归约入库

> 承接 #14 的 P4 外部加硬主线程。用户提供完整课堂视频 `20261021104500G214.mp4`
> （双曲线/椭圆教学，59.7 分钟，3 说话人）并确认"整堂课堂视频可以分析"。

#### feat: `experiments/r23_real_cognitive/`（课堂 S0 数据 + 采集协议）

- `thinkaloud-script.md`（NEW）：真人出声思维采集的完整操作手册——引导语 + 6 任务
  （A1-A3 酝酿 / B1-B3 反刍）+ 转写归约规则 + JSON 回传模板。S0 级真实认知数据源。
- `math-classroom-adaptation.md`（NEW）：数学课堂解题归约指南——候选类=解题途径、
  切段=换思路、权重=倾向伪计数、群体作答按学生分布估；含 M1-M6 建议题。S0 变体。
- `raw_transcript_lesson_G214.json`（NEW）：整堂课堂视频 ASR 完整转写（83KB，
  3 说话人，duration≈3584s）。
- `trajectory_human_G214.json` + `results_human_G214.json`：学生板演"求双曲线标准方程"
  段归约（定义法 2a=|PF1-PF2| vs 代数法 vs 焦点轴判断），判 **ρ_res>ρ_exp**（酝酿，
  osc_net=−0.92，4 步震荡）。
- `trajectory_human_001.json` + `results_human_001.json`：课堂截图（双曲线离心率）
  单条轨迹 M-dbl 归约，os_net≈−0.01 近中性。

#### 诚实边界（关键）
- 课堂是**讲授课**：教师主导、学生被点名引导式板演，**非独立出声思维**。
- 学生发言零碎（短句/配合式），转写噪声大（大量识别错误），归约权重为分析者推断。
- ⟹ 课堂板演归约**可靠性不足**，判 ρ_res>ρ_exp 仅作方法演示，**不能据此触发 X2 立项**。
- 与课堂截图 M-dbl 中性结论一致：课堂单线推导产生不了判别所需的多路候选权衡。
- 可靠 S0 采集仍走：① 单人独立出声思维（think-aloud 脚本）；② 多解/发散题独立求解出声。

#### 其他
- `formal/lake-manifest.json`：Lean lake 工程 manifest（v1.3.0，空 packages，formal 工程初始化）。

## v1.6.22 — 2026-10-24（本地未提交）

### 台账同步 + 新实验 pytest 覆盖（①+② 项）

**1. 七个新实验补上结构性测试（锁结构预言，不锁数字；仅标准库）**

| 测试 | 项数 | 锁什么 |
|---|---|---|
| `tests/test_a5_genericity.py` | 11 | genericity 的函数性 R=f∘Inv；峰序/截距随 sign(ρ−1) 分档且两档互斥；ρ_exp 共线反解；检查器有牙齿 |
| `tests/test_a7_core_min.py` | 11 | 传递闭包自反/传递；极小根集=可达∩根；孤儿记号被判不可达；η 登记为约定输入 |
| `tests/test_a13_first_case.py` | 9 | A_D 恒等式；U 非负；两叙事→互斥类别且峰序相反；D-CONV 方向与其 U 终值/峰高比自洽；四级 D-IDENT |
| `tests/test_b1_fixed_action_counting.py` | 7 | 四臂恢复四系数；闭合恒等式；β₀=κ_E·ρ_exp·(1−ρ) 三档异号；u₋/e₊ 记错破坏自洽 |
| `tests/test_b2_bridge_test.py` | 8 | 内点峰探测有牙齿；正常族误报<0.05、违反族≥0.9；**R31 是重言守卫、非 R43′**（边界随结果落盘） |
| `tests/test_b3_residual_inversion.py` | 10 | 推论 R43 的 ⟺ 双向（g 常数 ⟺ 残差无结构）；残差的仿射部分被吸收、**非线性部分不被吸收**；负对照无虚假结构 |
| `tests/test_r42_dscale_eta.py` | 11 | 判定表三支可达且完备；快振在粗窗口被平均、慢振不被平均；τ\*(η) 单调不增；读数与 η 无关 |

共 **67 项**（1.2 s 内跑完）；全套测试 = **368 passed**（含 numpy 环境；R29 单跑 78 s 需分开跑）。
`.github/workflows/xs-gates.yml` 的 pytest 步骤**自动覆盖**新测试，未改 workflow。

**2. 台账同步**：`CITATION.cff` 1.6.13 → **1.6.21**；README 顶部版本/测试数更新
（v1.6.6/199 → **v1.6.21/368**），并注明本机缺 numpy 时的 11 个 collection error 是
环境问题而非回归；`27-v3.4-boundary-map.md` **v0.3 → v0.4** + `boundary_map.json` v0.4。

**3. ⚠️ 新立礁石 R-043（τ\* 落点）—— 本次同步自己长出来的缺口**

R42 的扫描数据里已经写明「合理 η 区间 = 无」，但脚本的**结论段只说"判定稳定"**，
没有把这一事实写进结论。按诚实登记纪律补登记：

- 扫描中 τ\* 落在 (τ_fine, τ_coarse) 内的 η 区间**几乎为空**：A 同尺度与 C 结构性场景
  在 η 网格上**全区间失效**（τ\* 或 >τ_coarse、或退到最小窗口 2<τ_fine）；
  B 跨尺度弱场景仅 η=0.05 一格落在窗口内。
- **判定读数对 η 不敏感 ≠ η 的合理区间已定**。前者已验证，后者没有。
- 两条出路：① 加密 τ 扫描找回区间；② 若确实为空，τ\*/η 从**判定程序**降级为
  **解释性装置**（说明为何选这两个窗口），D-SCALE 证书须改写。
- 性质：纯内部可解，**不需要外部数据** ⟹ 目前清单上最便宜的一颗钉子。

**4. R-042 状态更新**：判定程序层已过（读数对 η 恒定、三类场景分落三支、g 不进判定），
但**不连带解除 B-B44**（η 取值原则）——该制约转由 R-043 承接。

**5. OPEN-4 · A/B 级噪声方差（纯仿真）**——补齐 16 §8.2/§13.2「有噪声未做」（`27-v3.4` §5 B7 行）：
合成 A_D 轨迹（复用 R30，θ_E=1.0，ρ∈{0.5,2.0}，4 调度族，N=200/配置）加 A 级（E 观测）/ B 级（U 观测）
噪声（σ_rel∈{0,0.001,0.01,0.05}，mode∈{A_only,B_only,both}），比较共线性法 vs 峰距法。脚本
`experiments/open4_noise_variance/`。
**关键发现（诚实，修正性）**：无噪声共线性法完胜（err~1e-15）；**A 级（E）噪声破坏共线性斜率**
（σ=0.05→err 0.91，差分放大）、峰序失守；**B 级（U）噪声影响小**（σ=0.05→err 0.02，峰序 0.99）。
⟹ 「共线性法最优」只在 E 观测噪声可忽略时成立。`27-v3.4` §5 B7 行 + 表 A 行 + 剩余开放格已更新
（OPEN-4 → 已闭合，剩 OPEN-3 / R-043）。

**6. G-26e · 非认知域对照（流行病防控）**——补 `26` §7.4 G-26e / `27` §5 B6 行：
A13 已在认知域（换工作决策）演示 D 系列判定程序；本项把**同一套程序机械移植**到**非流体、
非认知**域（流行病防控：E=累计已干预量，U=未决风险存量）。脚本 `experiments/g26e_cross_domain/`
（复用 A13 的 d_drill/d_conv/d_ident，未改判定程序）。
**结论（诚实限定）**：判定程序跨域可移植（演示层面）——流行病域上仍产出互斥类别
（钻研/空洞）、收敛方向、A/B/C/D 可辨识量；**非经验断言**，经验验证仍需多域真实数据。
登记：`26` §7.4 G-26e 行 + `27` §5 B6 行 → ✅ 已做。

**7. G-26f · P5–HMM 盲比行（预登记 + 文献）**——补 `26` §7.4 G-26f / `27` §5 B3 行：
把 `experiment_manual` §8.7 的 P5–HMM 链提升为盲比行，产出 `docs/xs-framework/28-v3.4-g26f-p5-hmm-blindtest.md`。
- **预登记**（查文献前）：P5 判"无特征尺度"（τ≈40 步穿越阈值）应对应 HMM 失效区间（P1）、且为硬失效（P2）、有特征尺度时 HMM 可用（P3）。
- **查文献**（HMM：E&Engquist 2005、Ming&Yue 2006、arXiv 2411/1803）：HMM 依赖尺度分离；中间尺度≈宏观网格 H → 有限误差；无尺度分离时 HMM/MsFEM 可能不收敛；LOD 型可在无尺度分离时仍收敛。
- **盲比结果**：**P1/P3 成立**（P5 作尺度分离探测器有效）；**P2 部分证伪**——"硬失效"过强，应修正为"精度受限/需降级或替代方法（LOD、细求解器）"，与 F2 优雅降级同型。
- **登记**：`26` §7.4 G-26f + `27` §5 B3 行 → ✅ 盲比通过；P5 判据使用说明补充（无特征尺度≠HMM 必然失效）；数值级验证待 G-26c。

**8. G-26a · ghost_gain_ratio 算子失配理论推导**——补 `26` §7.4 G-26a 行：
把 F2 的现象学注入参数 `ghost_gain_ratio`（默认 1.60）还原为界面算子结构量，产出
`docs/xs-framework/29-v3.4-g26a-ghost-gain-ratio-derivation.md`。
- **主结果**：`ghost_gain_ratio = 提升/限制算子（Restrict∘Lift）在副模态上的往返失配比 λ_i`。
  正常耦合（C2 过，λ_i=1）→ 失配项消失、守恒；F2（C2 失效，λ_i≠1）→ 注入净能量。
- **ghost_base_ratio** = 提升/限制算子离散误差（数值底噪）相对通量比率（patch test 残余，∝|gain|）。
- 与 R4 解析真值（r_i = −(g−1)·gain·δ_i − b·|gain|）、纯失配/纯底噪对照测试一致。
- **诚实边界**：理论推导非数值测量；λ_i 的数值待 G-26c 联合真实 MMS 求解器实测校准；
  内核当前是直接注入失配、未显式构造 Lift/Restrict 算子（其显式构造属后续工程）。

**9. GAPS X2 · CognitiveBridge.lean 立项评估**——产出 `formal/X2-cognitive-bridge-charter.md`，
更新 `formal/GAPS.md` X2 条目。
- **评估结论**：推荐**条件化立项**（两步）——① 先落实 ρ_res>ρ_exp 验证路径（B-B19 前置，
  补篇五 §1.6 的 R23：观察震荡期 U 走向判定），未触发前不写 Lean；② 触发后按验收标准写
  CognitiveBridge.lean（三命题以 ρ_res>ρ_exp 为显式资格条件，P1 取方向性版本，U≈0 数值阈值
  交仿真层）。
- **关键判断**：X2 是"待闭合的桥 + 待编码的命题"，不是单纯的"待编码命题"。跨层引用规则
  （未验证 ρ_res>ρ_exp 不得互引）要求先闭合 B-B19 前置，直接形式化会在签名层固化未经验证的桥。
- **风险登记**：B-B19 未闭合（高）、A_D 是建模假设 B-B22（中）、Ctx→Rat 数值映射换一致性论证
  路线（中）、定理 5.2′ 区间边界（中）、sign(0) 分支（低）。

**10. R23 验证路径 · ρ_res 与 ρ_exp 大小关系的可执行判定**——X2 步骤 1 落实：
产出 `docs/xs-framework/30-v3.4-r23-verify-path-rho-res-exp.md` +
`experiments/r23_verify_path/`（脚本 + results.json）+ `tests/test_r23_verify_path.py`。
- **操作协议**：① 令系统进入震荡（r=θ_E ⟹ Ψ≈0）；② 用补篇四程序 P̂ 读震荡期 U 走向；
  ③ U 降 ⟹ ρ_res>ρ_exp（许可跨层互引）、U 升 ⟹ ρ_res<ρ_exp、不变 ⟹ =。
- **数值验证**（固定种子，60 步 × 400 次）：三配置判别全部与真值一致；触底/触顶率精确复现
  补篇五 §1.5（I 100%/0%、II 0%/100%、III 0%/0%）。
- **稳健性**：配置 III（=）下单次 U 走向随机（53.5%/46.5%），判定须基于多次重复多数方向（>90% 阈值）。
- **诚实边界**：模型内自洽 ≠ 真实系统实证（A_D 是 B-B22 建模假设）；U 真实测量（程序 P̂ 对解释
  分布）是下游步骤，本脚本用模型真值；ρ_res/ρ_exp 单独标定仍不可行（B-B23，补篇九 定理 6.12）。
- **对 X2**：触发 CognitiveBridge.lean 立项的判据未达成（需真实系统实证 ρ_res>ρ_exp），步骤 2 不启动。

**11. 真实数据获取方案（B 线 P1）**——真实数据待办项的获取方案 + 首个候选实证：
产出 `docs/xs-framework/31-v3.4-real-data-acquisition-plan.md`（A 线 MMS / B 线认知，P0 定稿：
B 线先行、LLM 轨迹、R-043 排期）。
- **B 线 P1**（`experiments/r23_real_cognitive/`）：真实 LLM（Doubao）对 6 个构造歧义任务做真实
  多步推理，程序 P̂（补篇四 §3.2）+ 震荡期识别 + U 走向判定。
  **结果**：6/6 任务震荡期 U 下降 → 判 **ρ_res>ρ_exp**（酝酿型），X2 立项判据**初步触发**。
- **方法修正**：震荡期识别初版误加"U 平稳"条件，把酝酿/反刍误筛为无震荡（只能判出空转）；
  修正为"多候选并存（r̂≥2）即震荡，U 走向为判别量"（符合推论 5.3）。
- **诚实边界**：构造偏倚（全收敛型样本、无反刍对照）、单模型单次采样、自生成轨迹（S1 硬度）、
  不构成一般结论。正式立项需 **P4 加硬**（多模型、真人 think-aloud/编辑历史、含反刍对照）。
- **R-043**：确认纯内部可解，排期（不占真实数据线）。

**12. R-043 · τ* 落点（η 合理取值区间）**——纯内部可解，已解决：
产出 `experiments/r43_tau_star/`（r43_tau_star.py + results.json + README），
登记 reef_shturl.json R-043 + boundary_map.json B3 格。
- **结果**：加密 η 网格（200 点）后，三类场景在默认窗口 (8,128) 下均有非空 τ* 落带 η 区间
  （A [0.031,0.046]、B [0.031,0.071]、C [0.061,0.086]）⟹ **R-042 的『几乎为空』是网格假象**
  （0.02→0.05 跳变错过落带区间）。
- **B-B44 降级**：从"η 取值原则开放"降为「已给出各场景经验区间、待锚定取值」。
- **诚实边界**：η 合理取值**场景相关**（三类交集为空：A 上界 0.046 < C 下界 0.061），区间窄，
  仅本仿真三类场景/固定种子下成立；D-SCALE 判定读数仍与 η 无关（R-042 结论不变）。

**13. P4 B 线加硬 · 反刍对照（双向判别验证）**——补齐 P1 缺失的反刍方向：
产出 `experiments/r23_real_cognitive/`（trajectory_rumination.json + results_rumination.json +
P4-hardening-protocol.md，README 整合）。
- **结果**：4 个反刍诱导任务（理发师悖论/无信息选门/自指/矛盾数）全部震荡期 U 上升，
  判 **ρ_res<ρ_exp**（反刍型）。与 P1 的 6/6 ρ_res>ρ_exp（酝酿型）形成**干净双向对照**：
  正常可收敛任务→酝酿（ρ_res>ρ_exp），悖论/无信息/矛盾任务→反刍（ρ_res<ρ_exp）。
- **意义**：验证 R23 判别在两个方向都有效，且判别与任务结构一致（推论 5.3 语义成立）。
- **诚实边界**：仍为单模型（Doubao）自生成轨迹（S1）、样本小（P1 6 + P4 4 = 10 任务）；
  正式触发 X2 立项仍需 **P4 外部加硬**（多模型轨迹、真人 think-aloud/编辑历史），
  采集协议见 P4-hardening-protocol.md（需用户外部采集回传）。

**14. Doubao 扩展样本 + X2 触发评估**——把 B 线扩至 18 任务：
产出 `experiments/r23_real_cognitive/`（trajectory_extended.json 新增 8 任务 +
trajectory_all.json 合并 18 + results_extended.json + results_all.json，README 更新）。
- **汇总（18 任务）**：ρ_res>ρ_exp 12/18（67%）、ρ_res<ρ_exp 6/18（33%），多数酝酿。
- **关键发现**：判别捕捉实际 U 走向、**不依赖任务标签**——悖论任务 R7/R8/R9 判酝酿、
  正常权衡 T8 判反刍 ⟹ 判别器独立性增强，但"悖论→ρ_res<ρ_exp"简单叙事被修正。
- **X2 触发评估**：酝酿方向仍多数成立，触发判据保持**"初步满足但不充分"**；反刍对照
  存在但不稳定。正式立项仍需 P4 外部加硬（多模型/真人）。登记 X2 charter 步骤 1 更新。

## v1.6.21 — 2026-10-24

### B2 · 桥检验（R31：A_D 一致性守卫）——方法 + 脚本准备（不声称经验验证）

> 对应 `18-v3.3-formal-framework-charter.md` §6.2 **B2**（缺口 B-B19 / R43′）。
> B 组需外部数据；本项做「数据到达前的可执行准备」。

- **范围澄清（诚实）**：R43′（β₀=−ρ_exp·Ψ(k*_U)，独立标定 Ψ）是真桥检验，但需补篇四
  独立标定程序 + 真实认知数据，**不在合成数据方法准备范围内**。本项实现其**重言形式 R31**
  （β₀=−ρ_exp·ΔE(k*_U)）作为 **A_D 一致性守卫**——验证数据是否满足 A_D + 内点峰自洽。
- **方法**：合成 A_D 轨迹（钻研型，U 内点峰）→ 回归 `(ρ_exp,β₀)` → 峰点检验
  `|ρ_exp·ΔE_c+β₀·dt|/scale`。两族对比：正常族（纯 A_D）vs 违反族（A_D⁺ 边界层注入
  `ΔU+=κ_E·ρ_exp·(1−ρ·g(U))·dt`，g(U)=g_strength·U）。
- **脚本**：`experiments/b2_bridge_test/b2_bridge_test.py`（仅标准库，确定性）。
- **结果**（`b2_results.json`）：**正常族误报率 0.000**（rel_viol≈1e-4）；**违反族误报率 1.000**
  （g∈{0.4,0.8,1.2}，rel_viol 0.95–1.0）。R31 判别能力验证通过。
- **量纲修正（本项暴露并修复）**：峰点检验初版用 `β₀·2·dt`（应为每步 `β₀·dt`），致正常族误报 1.0；
  修正后正常族 0.000。违反族由"ε·U 状态依赖"改为"A_D⁺ 边界层注入"（理论自然、稳定 →1.0）。
- **诚实边界**：① R31 是重言形式守卫，通过≠桥 B-B19 成立；② R43′ 真桥检验未做（需独立 Ψ）；
  ③ 合成数据非经验验证；④ 需 U 有内点峰。
- **登记**：`18-v3.3` §6.2 B2 行标注「R31 守卫方法+脚本准备已完成 v1.6.21，R43′ 真检验仍属 B 组」。

## v1.6.20 — 2026-10-24

### B1 · 固定动作计数注入（R33 / 命题 6.3）——方法 + 脚本准备（不声称经验验证）

> 对应 `18-v3.3-formal-framework-charter.md` §6.2 **B1**（缺口 B-B45 / R45，闭合 C 级）。
> B 组需外部数据；本项做「数据到达前的可执行准备」——实现并验证四臂闭合分析流程。

- **方法（命题 6.3）**：臂1（a=1,b=0）测 `(e₊,u₊)`，臂2（a=0,b=1）测 `(e₋,u₋)` ⟹ 一次闭合
  `θ_E=e₋/e₊`、`θ_U=u₋/u₊`、`ρ_exp=u₊/e₊`、`ρ_res=u₋/e₋`、`κ_E=e₋`、`β₀=ρ_exp·e₋−u₋`。
- **脚本**：`experiments/b1_fixed_action_counting/b1_fixed_action_counting.py`（仅标准库，确定性合成）。
- **结果**（`b1_results.json`）：3 配置（drill/converge/boundary，借 R27 CONFIGS）四臂闭合 vs
  轨迹回归自洽 **rel_ρ≈1e-14、abs_β≈4e-14**，全部 PASS。
- **判据修正（本项发现）**：`β₀` 真值可为 0（如 boundary 配置），**相对差会因除零失效**，
  改用**绝对差**判据 `abs_β < 1e-6`。
- **诚实边界**：① 合成数据非真人实验（四臂真实测量含噪声，加噪稳健性留作扩展）；
  ② `κ_E=e_minus` 为仓库惯例（a13 make_cfg），理论定义若有出入需核对；
  ③ C 级闭合（R45）真正达成需真人实验数据（B 组性质）；④ 回归依赖 A_D 成立。
- **登记**：`18-v3.3` §6.2 B1 行标注「方法+脚本准备已完成 v1.6.20，真人测量仍属 B 组」。

## v1.6.19 — 2026-10-24

### B3 · 残差反演 g（R32）——方法 + 脚本准备（不声称经验验证）

> 对应 `18-v3.3-formal-framework-charter.md` §6.2 **B3**（缺口 B-B43 / R32）。B 组需外部数据；
> 本项做「数据到达前的可执行准备」——实现并验证残差反演 `g` 的方法本身。

- **方法**：A_D⁺（15-v3.2 §5.7 推论 R43）下 `R(k)=U−ρ_exp·E=α+βk−κ_E·ρ_exp·ρ·∫₀ᵏ(g−ḡ)dt`；
  对残差**去仿射趋势 → 中心差分 → 除以 `κ_E·ρ_exp·ρ`** 即反演 `g(U)`（差常数平移，规范自由度 B-B23）。
- **脚本**：`experiments/b3_residual_inversion/b3_residual_inversion.py`（仅标准库，确定性合成）。
- **结果**（`b3_results.json`）：4 种单调增 g(0)=0 的 g（linear/square/tanh/saturated）反演
  **corr=1.000000、norm_rmse≲1e-3、scale≈1**；负对照（g≡0、g≡常数）反演 std≈1e-14（**无虚假结构**）。
- **诚实边界**：① 方法验证非经验验证（"真实残差是 g 信号"不由此声称，属 B4 跨人群）；
  ② 反演恢复 g 的**变异/形状**，绝对水平不可辨识（B-B23 常数平移）；③ 振幅比例依赖
  `κ_E·ρ_exp·ρ` 独立标定（关联 B2/B-B19）；④ U 动态范围较窄（受 ρ·g 平衡约束），更宽范围
  与加噪稳健性留待后续扩展。
- **登记**：`18-v3.3` §6.2 B3 行标注「方法+脚本准备已完成 v1.6.19，套用真实数据仍属 B 组」。

## v1.6.18 — 2026-10-24

### A9 · Lean4 形式化接入（candidates/DualConstraint.lean → Formal/）

> 用户回到有仓库的原电脑并装好 Lean 工具链后，执行 A9（Lean4 全内核形式化）的核心接入步骤。
> 对应 `formal/GAPS.md` 条目 C1（认知双约束 / R18 形式化）的接入步骤。

#### 环境
- 本机装好 elan 4.2.4，`lean`/`lake` 可用。
- **工具链（最终决策）**：`lean-toolchain` 原指向 `leanprover/lean4:v4.14.0`，但本机从
  `releases.lean-lang.org` 下载 v4.14.0 反复因 **SSL 连接被重置** 失败；stable（v4.35.0）可下载，
  故先用 `elan override` 以 v4.35.0 编译验证（Lean 4 基础语法 v4.14→v4.35 兼容，`lake build` 全绿）。
  **用户 2026-10-24 拍板：`lean-toolchain` 升级到 `leanprover/lean4:v4.35.0`**，移除 override 后
  由 lean-toolchain 直接生效，`lake build` 复核仍 8 jobs 全绿。今后项目统一以 v4.35.0 构建。

#### feat: `formal/Formal/DualConstraint.lean`（接入）
- 将 `formal/candidates/DualConstraint.lean`（认知双约束 / R18 的 Horn skeleton 编码，
  原"候选，未编译"）**单独编译通过**后移入 `Formal/`。
- `Formal.lean` 增加 `import Formal.DualConstraint`。
- `lake build` 全项目（Basic/Axioms/Corollaries/AttackSurface/DualConstraint/Formal，8 jobs）**编译通过，无 sorry**。
- 编码内容：两类失效（外/内）的可分离性与判别签名 R18-a（DC3–DC8 公理 + T1–T6 可导出定理）；
  定量签名（数值不等式）与依赖 B-B19 的桥按 `GAPS.md` 登记为 X1/X2 不编码。

#### docs
- `formal/GAPS.md`：C1 条目状态由「候选，未编译验证」改为「已接入」。
- `formal/candidates/DualConstraint.lean`：按接入步骤"移入"，从 candidates 移除（避免重复），
  内容见 `Formal/DualConstraint.lean`。

## v1.6.17 — 2026-10-24

### A13 · 第一例对应探索（换工作决策）

> 用户指令：「A5 完成后接着 A7、A13」。本轮完成 A7（上一条目）与 A13（本条目）。
> 对应 `18-v3.3-formal-framework-charter.md` §6.1 触发点 A13——**投影空间与现实空间的第一个正式接触点**。

#### feat: `experiments/a13_first_case/`（a13_first_case.py + a13_contrast.json + README.md）

对"换工作决策"单一认知案例，用 D-DRILL / D-CONV / D-IDENT 跑完整判定流程，输出
**框架判定 vs 无框架人工判断**对照记录：
- 场景 A（越查越笃定，ρ=2）：D-DRILL **钻研相**（U 峰在前），D-CONV **收敛方向**（U 降到峰高 7%）
- 场景 B（越查越乱，ρ=0.6）：D-DRILL **空洞相**（E 峰在前），D-CONV **发散方向**（U 保持峰高 96%）

**对照点**：无框架把"清楚/乱"当程度副词并常归因于个体特质；框架把它判读为互斥的钻研/空洞相位 +
可定位的可辨识量（ρ_exp 斜率 A 级、sign(1−ρ) 截距 B 级、κ_E 需外部单位 C 级、四单位系数 D 级伪目标）。
可作论文示例 / Graphical Abstract 素材 / 审稿演示案例。

**诚实边界**：演示性映射，非现实测量（真实决策无客观可测 (E,U)，叙事→轨迹是构造）；
D-CONV 只判收敛/发散方向，不裁真伪（真/伪收敛需外部证伪信息，超出本演示）；跨现实验证属 B 组 B4。

#### docs

- `18-v3.3-formal-framework-charter.md`：A13 触发点→已完成。

## v1.6.16 — 2026-10-24

### A7 · 核心最小化（47 记号 → 极小核心 + 导出表）

> 用户指令：「A5 完成后接着 A7、A13」。本轮完成 A7（核心最小化）。
> 对应 `18-v3.3-formal-framework-charter.md` 待办 A7。

#### feat: `experiments/a7_core_min/`（a7_core_min.py + core_table.json + README.md）

对 `symbol_registry.json`（47 记号）做依赖图机检压缩：
- **核心 = 11**（根 10：A_D/C_in/F/M/O/O_all/Ω/S/X/φ + 约定输入 1：η）
- **导出 = 36**，全部可由核心沿 `derives_from` 到达，**展开校验通过（不可达 0）**
- **核心压缩率 76.6%**；L2 的 15 个动力学参数全部由 A_D 导出（体现 A_D 唯一承重假设地位）

**诚实边界**：压缩率是**依赖结构层面**，非**概念独立性层面**（10 根是否彼此独立是概念判断，
不在图论机检范围）；未改动任何记号定义与层归属，未触碰冻结文本。

#### docs

- `18-v3.3-formal-framework-charter.md`：A7 状态→已完成。

## v1.6.15 — 2026-10-24

### A5 · generics 一般性证明（代数证明 + 原型枚举）

> 用户指令：「继续下一个待办」⟹ 执行 A5（generics 一般性证明）。
> 对应 `18-v3.3-formal-framework-charter.md` §1.4（genericity 等价式）、判据 C3.1、待办 A5。

#### feat: `experiments/a5_generics/`（a5_genericity.py + a5_results.json + proof.md）

在 A_D 原型上证明并机器核验四条分类规则的 genericity（输出只依赖不变量）：
- 规则 1 楔形分类（D-DRILL）：Inv=sign(ρ−1)，R∈{drill,hollow,none}
- 规则 2 峰序（P3）：Inv=sign(ρ−1)，R=U峰前后
- 规则 3 四级可辨识性（D-IDENT）：Inv=(ρ_exp, sign(ρ−1))，自足两级
- 规则 4 共线性反解 ρ_exp：Inv=ρ_exp，R=OLS 斜率

**机器验证**：720 实例（4 调度 × 3 θ_E × 3 ρ_exp × 5 ρ × 2 E0 × 2 U0），
四条规则同不变量组内输出恒同，**零违反**；共线反解跨实例 max rel err 1.48e-16。

**形式精确化（重要）**：章程 §1.4 的 ⟺ 双向式对细粒度不变量过强
（楔形分类中 ρ=1/3 与 ρ=0.5 都判 hollow）。正确的 genericity 形式是**函数性** R=f∘Inv
（f 不必单射）。证明文档 `proof.md` 详述。

**诚实边界**：代数 + 原型枚举，**非 Lean 形式化**（Lean 归 A9，本机无工具链）、非现实断言。

#### docs

- `18-v3.3-formal-framework-charter.md`：A5 状态→已完成、判据 C3.1→原型内已证、
  可说/不能说清单第 4 条加 v3.4.1 收尾登记（限定"已在 A_D 原型内证成"）。



### A10 收尾（C1.7/C1.8 接入 CI 强制门）+ R-042 D-SCALE 参数敏感性扫描

> 用户指令：「直接动手执行」A 组待办。本轮完成第一优先 A10 收尾（接 CI 强制门）与第二优先 A12/R-042（D-SCALE 扫描）。
> 对应 `18-v3.3-formal-framework-charter.md` §6.1（A10 状态说明、优先级表）与 `reef_shturl.json`/`distinction_registry.json`。

#### feat: `.github/workflows/xs-gates.yml`（A10 · CI 强制门）

把 XS 合法性机检接成 CI 合并门槛，三道强制门任一 error 即阻断：
`xs_trace_gate.py`（C1.7 推导链可达性 + C1.8 冻结窗哈希比对）、`xs_formal_checks.py`（C1/C2 记号与区分登记一致性）、`pytest tests/`（全套数值原型）。

**本地等价验证全绿**：两个校验器 error=0（唯一警告为已登记的 `η` 无根断言）；`pytest` **301 passed**；
坏样本注入验证强制门确实拦截（G001 被捕获）。
**诚实边界**：CI 门需推送到 GitHub 并启用 Actions 后才是真正的合并门槛；本地仓库无 GitHub Actions 运行时，以上为本地等价验证。

#### feat: `experiments/r42_dscale_eta/`（R-042 · D-SCALE 对 η 的参数敏感性）+ 登记

- `r42_dscale_eta_sensitivity.py`：三类场景（同尺度/跨尺度弱扰动/结构性扰动）× η∈[0.005,1.0] 网格，
  固定观测窗口 τ_fine=8/τ_coarse=128、无噪声（框架判定为无噪声语义），固定种子 20261024，逐位可复现。
- 结果：三类场景判定读数**正确区分且对 η 稳定**——A→same_scale_coupling（R 两尺度≈0）、
  B→cross_scale_weak（R_fine=0.195>0, R_coarse≈0）、C→structural（两尺度>0）；
  g 不进入判定流程 ⟹ 支持 §4.4.1 B-B44「弱」拆分。
- **登记**：`reef_shturl.json` R-042 状态→「已扫描」；`distinction_registry.json` D-SCALE status→**established**（附 B-B43/B-B44 残余说明）。
- **方法学澄清（诚实登记）**：τ_fine/τ_coarse 是固定观测窗口，τ* 为应落在两者间的分离尺度，**判定读数与 η 无关**；
  η 决定 τ* 是否落在窗口之间（有效区间约束）——不是"不稳定"，需在证书写明的边界。

#### docs

- `18-v3.3-formal-framework-charter.md`：A10 状态收尾（v3.4.1 登记）、优先级表 A10/A11/A12 三行改完成、
  可说/不能说清单第 6/7/10/12 条加 v3.4.1 收尾登记。



### R33 C 级闭合实测（R45 核验）+ R34 ρ=1 顶点容差带（OPEN-1 的解）

> 用户指令：「解掉 OPEN-1（容差带）与 R45（C 级实测核验）」。
> 对应 `27-v3.4-boundary-map.md`（边界地图 v0.1）的两个空白格。

#### feat: `experiments/r33_clevel_closure/` + `tests/test_r33_clevel_closure.py`（6 项）

命题 6.3「固定动作计数闭合 C 级」的实测核验（R45 判据）：四臂
（固定展开 1 次 / 固定消解 1 次 / 自由调度序列 / 回归 `U = α + ρ_exp·E + β·k`）
在 A_D 动作计数原型上对账。**36 行（3 θ_E × 3 ρ × 4 调度族）T1–T4 全过**，
残差为浮点噪声；**负对照有牙齿**（b 记错破坏 T2、不动 T1；ρ=1 退化点豁免已登记）。
⟹ R45 原型内成立，B-B42/B-B45 由「未做」改「已核验（原型级）」。

#### feat: `experiments/r34_vertex_band/` + `tests/test_r34_vertex_band.py`（6 项）

OPEN-1（B 级符号判读在 ρ=1 顶点退化）的解：**B 级判读禁区（顶点带）:=
{|β̂| ≤ z·SE(β̂)}**，β̂/SE 来自 `U = α + γE + βk` 回归，z 约定性（默认 1.96）。
判读升级为**三值**（钻研方向 / 空洞方向 / 顶点带），不强行二分。
结构检验：无噪声禁区恰为 ρ=1；有噪声带宽随 σ 增宽、留空集中顶点邻域；
**G_λ 下 (E,U) 逐点不变、判读不变（U5 通过）**。

#### docs

- `27-v3.4-boundary-map.md`：OPEN-1 关闭、R45 核验通过、B 行例外条款改写、
  §5 交叉索引与总评更新；剩余开放格 OPEN-3（开闭二分）/OPEN-4。
  **同日升 v0.2（拆表）→ v0.3（约定化）**：表 A 分级结构【定理】（四级不随实例变化）
  与表 B 可读性条件分开；C 行 = 辨识层【定理】+【原型确认】（R33）+ 自足层【约定】；
  **OPEN-2 从【开放】（等待判据）改登记为【约定】（等待声明）**——约定 A（观测工具
  算界内）下 C 级自足为定义上的真，约定 B 下为相对于观测协议的条件性结论；
  框架不裁定约定、只要求显式声明（与 η、z 并列的约定性输入）；R45 与 OPEN-2
  登记为同格两面；D-IDENT 证书附款按约定声明生效。
- `boundary_map.json`：open_cells / gap_index / r041_table 同步。
- README 缺口表 B-B42/B-B45 改「已核验」；与代码层对应新增两行。

---

## v1.6.12 — 2026-10-23

### XS 章程 v3.4 修订：前置生成理论 P 的悬空登记、三类探索谱系、A 组优先级重排

> 用户指令：「修正并细化原 D2 界面规则，不再是单纯一刀切的出界切割；把**悬空但起作用的
> 前置生成理论 P** 显式登记、隔离、约束，同时重排 A 组工作优先级，配套
> "边界地图优先"的完备性判据。」

#### docs: `18-v3.3-formal-framework-charter.md`（**v3.4 修订版**，五处改动）

**① D2 界面规则改写（§1.3）**：旧 D2 = **切断**来源讨论；新 D2 = 来源假说 P
**单独悬空登记、不进内核**，但保留一条**受控、可操作**的通道去探索两边的对应关系。
旧文保留（U1），加取代标注。

**② 新增 §2.5「前置生成理论 P 与判别内核的二分」**：
- **P1** P 未被认定：不进 A1–A11，**不参与 C1/C2/C3**，框架自洽性不依赖 P 为真；
- **P2** P 的功能是**生成性**（为什么有权构造投影空间、为什么有这些变量）；
- **P3** 内核的功能是**判别性**（只在给定投影空间的前提下定义区分与读法约束）。
- 三条推论：即使 P 假，内核仍可自洽（类比欧氏几何不依赖物理空间平直）；
  但若 P 无现实锚点，内核退化为纯记号游戏；投影↔现实的关系是**开放研究问题**。
- 含**薛定谔类比**（德布罗意假说 = P，能级 = 判定输出，光谱比对 = 对应/约束探索）
  及其**辖域声明**：类比只说明结构有先例，**不构成对 XS 的辩护**（XS 一个接触点都没有）。

**③ §6.1 A 组按依赖顺序重排**（原则：**先拿到可信、完整的边界地图**）：
- **够用判据**：边界地图完成 + C1.7/C1.8 机检落地 = 可以启动 A13 对应探索，**不必等 generics 与 Lean**；
- 第一优先：A10（C1.7/C1.8 机检）+ A11（边界地图，R-041）；
- 第二优先：A12（D-SCALE 参数敏感性，R-042）；
- 第三优先：A5（generics）、A7（核心最小化）；
- **可推迟**：A9（Lean4）——边界测绘若改了规则，形式化代码会反复重构；
- **触发点 A13**：第一例对应探索（换工作决策案例，D-DRILL/D-CONV，输出框架判定 vs 无框架人工判断）。

**④ 新增 §6.4 三类探索谱系**（新 D2 那条受控通道的具体内容，强度递增）：
**对应探索（弱）**——区分判定在案例上是否稳定有差异，**不验证 P**；
**边界探索（中）**——可读/不可读边界是否与现实认知的限制结构对齐；
**约束探索（强）**——形式不可能性能否传递为现实不可能性（示例：`(κ_E,ρ)` 不可分离）。
每项写明**目标**与**不验证什么**（防"受控通道"退化成口袋）。

**⑤ §4.4.1 D-SCALE 两条制约的强度拆分**：
**B-B44（η）判为弱**——`τ*` 存在性只依赖 η 落在合理区间；
**B-B43（g）判为弱**——判定只需粗细两尺度残差对比，**g 不进入判定流程**（它阻塞的是 R32 反演，不是 D-SCALE）。
⟹ 此前把两者捆为同一强度，是**归类错误**。但**拆分是论证不是检验**，状态仍为 provisional（需 R-042 扫描）。

**§7 代价清单新增 2 条**：P 悬空的代价（可能根本没有验证）、受控通道的滑坡风险。
**§9 不能声称的新增 5 条**（第 8–12 条），含"不能说 C1.7/C1.8 已落地"。

**§6.3 边界裁定存疑（待裁决，不擅自改判）**：新 D2 下 **A_D 的经验地位**恰好落在边界上——
说它出界（需本体承诺，属 P 辖区）与说它在界内（与边界/约束探索性质相同）都能论证。
本稿**仍按出界处理**，理由是需要本体承诺的部分占主导；若改判则成为 B5 并与 §6.4.3 合并。

#### docs: `reef_shturl.json`（NEW · 礁石台账）

与 B-Bxx 缺口表分工：缺口表记"框架的缺口"，礁石台账记"绘制边界地图过程中的待解问题"。
- **R-041**：可辨识分级是否随开闭实例变化（D/C/B 三级逐条核查清单 + 预期，产出二维边界表）；
- **R-042**：D-SCALE 对 η、g 的参数敏感性（含强度拆分、扫描程序、两个分支结果的处置）。

#### docs: `distinction_registry.json`

新增两字段：`boundary_dependence`（对应 R-041）、`parameter_sensitivity`（对应 R-042）。
**当前六条区分全部"待测绘/未测"**，校验器报 **D010 warn** 作为待办的**显式回声**。

#### src: `src/xs_trace_gate.py` + `tests/test_xs_trace_gate.py`（**20 项测试**，仅标准库）

把此前**只有设计、没有代码**的 C1.7 / C1.8 变成机检项：
- **C1.7 推导链溯源**：每个非根记号必须沿 `derives_from` 可达根记号。
  **根 = 公理级记号 ∪ 已声明假设**（采自 C1.7 原文"回溯到公理**与已声明假设**"）。
  违规码 G001（未登记）/ G004（白名单过期）/ G005（白名单条目缺 reason 或 todo）/ G000（无根记号）。
- **C1.8 冻结窗**：`frozen_manifest.json` 登记 11 份 v3.1 主文+补篇+审计文本的 sha256，
  逐条比对。违规码 G002（被篡改）/ G003（被删）/ G006（清单缺失→只是口头声明）/ G007。
- **哈希按"换行归一化为 LF"后计算**（`sha256_of(normalize_newlines=True)`）。
  原因：git 在 Windows 上按 `core.autocrlf` 检出时会把 LF 写成 CRLF，
  同一份文本在不同平台字节流不同 ⟹ 直接哈希会**全库假阳性**。
  v3.4 首次生成清单后的干净克隆验证**实测报出 9 条 G002**，全部是换行差异。
  归一化后**文字改动仍会被抓到**（已加两条测试分别锁住"忽略换行"与"仍抓内容改动"）。
- 已并入主校验流水线：`python src/xs_formal_checks.py` 会一并运行。

#### docs: `trace_gaps.json`（NEW）· `frozen_manifest.json`（NEW）

**首次全库回归扫描结果**：47 条记号中 38 条非根，**24 条回溯不到根**。
逐条核对 `definition` 原文后：
- **10 条是登记表漏写 `derives_from`**（`a,b,e₊,e₋,u₊,u₋,g` → `A_D`；
  `P_flx` → `P_lvl`；`Q_flx` → `Q`；`E` → `Psi`）⟹ 已补；
- **13 条因连锁恢复可达**；
- **真正无根的只剩 1 条**：`eta`（约定性输入，B-B44）⟹ 登记在 `trace_gaps.json`。

#### ⚠ 本轮第二处自我更正（诚实登记）

初稿把 `E`、`P_flx`、`Q_flx` 也判为无根，理由写作"通量不由存量推导"。
**核对定义原文后该理由不成立**：`P_flx` 的定义就是"P 的通量形式"，`E` 的定义就是
`E(k+1)−E(k)=Ψ(k)`——**错在把"标定通道独立"当成了"定义独立"**。
校验器的 **G004（白名单过期）当场抓出了这三条**，已从白名单移除。
⟹ 这是校验器**第二次**抓出我自己的错误，且抓的是**登记行为本身**
（往白名单里塞了一条不该塞的）——白名单因此没有退化成"让洞消失的口袋"。

#### 与方案原稿的一处出入（如实记录）

方案原稿称「`derivation_trace.py` / `frozen_gate.py` 代码就绪，尚未并入流水线」。
**核对仓库后：这两个模块不存在**（`src/` 下只有 `xs_formal_checks.py` 及早期实验脚本）。
v3.4 **新建**了 `src/xs_trace_gate.py`。

#### 不能声称的

C1.7 / C1.8 **未完全落地**：全库回归扫描的**结果登记**未完成、**CI 强制门未接**
⟹ "所有结论可回溯至公理、冻结基准无静默篡改"仍是**目标**而非已核验的框架性质；
P 未被检验；D-SCALE 仍为 provisional；R-041 / R-042 未做；A_D 的经验地位边界裁定待裁决。

---

## v1.6.11 — 2026-10-23

### XS 形式框架章程（v3.3 纲领）：定位变更 + 三条合法性判据的机器化

> 用户指令：「假如现在我构建的不是一套关于认知的理论，而是一个用于组织认知过程描述的形式框架，
> 我们仅仅构建 XSM，边界不外溢。它的合法性不来自经验验证，而来自内部一致性、术语明确性和
> 跨实例的可复用性。它的意义取决于是否有人能用它做出在框架外部无法做出的区分……
> 现在我们可以做什么？我们已经做了什么？」

#### docs: `docs/xs-framework/18-v3.3-formal-framework-charter.md`（NEW）

定位变更：XS **不是**关于认知的理论，而是**组织认知过程描述的形式框架**。

| 维度 | 旧定位 | 新定位 |
|---|---|---|
| 合法性来源 | 经验验证 | **内部一致性 + 术语明确性 + 跨实例可复用性** |
| 核心产出 | 关于认知的断言 | **可复用的区分** |
| 缺口性质 | 待现实证据闭合 | 多为**声明不做的事**或**待补的判别力证书** |

**新增 D2 · 界面规则（边界不外溢）**：框架只对 XSM 实例说话；
凡涉及「实例怎么来的」的问题一律出界 ⟹ **B-B22 从「最高价值缺口」改判为「声明不做的事」**。

**XSM 实例的形式定义**（§2）：闭实例（`s(k)` 可观测，C 级可闭合）vs 开实例（不可观测，C 级不可分离）。
⟹ 「为什么有些量读不出来」从**技术故障**变成**实例类型的性质**。

**三条判据的可操作化**（§3）：C1.1–C1.6、C2.1–C2.8 全部机检；C1.7 / C1.8 / C3 列为未完成。

**判别力四件套**（§4.1）：非平凡实例 / 可执行判定程序 / **外部必要性** / 状态。
**六条区分**（§4.2）：D-IDENT / D-SCALE / D-CONV / D-DRILL / D-LOCAL / D-FAIL。
§4.3 用「为什么随便挑两个量不行」示范外部必要性（D 级是显式反例）。

**下一步**（§6）：A 组 10 项**不需要外部证据**（A1–A4 本轮已做；A5 generics 一般性证明 /
A6 判别力证书 / A7 核心最小化 / A8 术语词典机器化 / A9 Lean / A10 C1.7/C1.8 机检化）；
B 组 4 项**需要外部证据**（B1–B4，**不阻塞** A 组）；C 组 3 项**显式出界**。

**§7 代价清单**（诚实）：对现实的解释力归零；意义完全押在判别力上；
判据三有滑坡风险（"在自己造的类里自洽"），防线是 B4 与 §4.1 件③。

#### docs: `docs/xs-framework/symbol_registry.json`（NEW，47 条记号）

每条记号的：层、种类、量纲（含"无量纲"）、定义、依赖、**可辨识级别**（A/B/C/D）、
**允许的读法**、**禁止的读法**、定义出处。
外加 `legacy_conflicts` 2 条：**C1「Ψ / U 符号冲突」**、**C2「Q / Q_prune 术语二义」**。

#### docs: `docs/xs-framework/distinction_registry.json`（NEW，6 条区分）

每条区分的四件套。用户点名的两条均已登记：D-IDENT（读得出 vs 读不出）、
D-SCALE（同一尺度下的耦合 vs 跨尺度的弱扰动，provisional）。

#### src: `src/xs_formal_checks.py` + `tests/test_xs_formal_checks.py`（24 项测试，仅标准库）

把三条判据变成可机检项：依赖无环（S013）、无自指（S007）、层归属唯一（S002/S006）、
记号唯一（S012）、量纲显式（S004）、可辨识级别必填（S008）、允许/禁止读法不相交（S011）、
遗留冲突有处置（S014/S015）；区分侧：非平凡实例 ≥2 且互异（D003）、外部必要性非空（D005）、
不可执行不得标 established（D004）等。

**校验器首次运行就抓出本登记表自己的 17 处违规**（量纲缺失、派生缺失、未登记引用、区分侧拼写错误）。

#### ⚠ 本轮一处自我更正（诚实登记）

v3.2 初稿把 `Ψ`/`U` 的记号冲突写成**「新发现的第九个病根 S9」**——实为 **v3.1 补篇四 §2.4 已登记**。
**作者在自己的新文章里重犯了一个已经修过的错。**
⟹ 这正是新定位的实证：**靠人读，记号冲突会漏检**；所以需要登记表与校验器。
已更正 `15` / `16` / `17` / README / CHANGELOG；并在补篇三 §3.3、补篇五 §2.2 加了
取代/交叉引用标注（U1，原文保留）。

#### 不能声称的

XS 框架不描述真实认知过程；合法性来自三条形式判据而非现实检验；判别力是**论证**不是使用记录；
generics 未证（A5）；跨现实实例未做（C3.4）；D-SCALE 标 provisional；校验器不覆盖 C1.7 / C1.8。

---

## v1.6.10 — 2026-10-23

### XS 框架 v3.2 通俗版：一份能直接转给非专业读者的文章

> 用户指令：「我需要一份象我提供的文章一样的文字，普通人也能看懂一些。」
> 参照文：《我对认知的个人理解》（第一人称、口语化、每节配真实案例、
> 明确标注「不是客观真理」、含翻车记录与边界清单）。

#### docs: `docs/xs-framework/17-v3.2-plain-language-essay.md`（NEW）

**内容与 `16` 严格一致，只换表达。** 十节：

1. 认知不是"想得对"，是"想的过程有形状"（含四个作用与四个性质）
2. 先拆别人的框架：心理学 / 哲学 / 脑科学 / 认知科学 / 玄学 **各有边界，皆非万能**
3. **这套框架的四次改写**——每段配真实案例：
   胶带机器（亲手废掉唯一证伪路径，误报率 → 100%）→ 开始归因（21 项里 13 项是自己挖的坑）
   → 推翻自己（标定式是循环的 / 桥被写成恒等式 / 证否途径其实早就有效）
   → 现在在哪
4. **一条公式的大白话版**（库房比喻）与"为什么一条公式就够了"
5. 五种读法（哪两种稳、哪两种脆；峰距法 36%–76.5% 偏差的真正原因：越层读取）
6. **形态**：楔形图（钻研与空转是同一楔的上下两半，互斥）
7. **能读什么不能读什么**（"有些数天生读不出来"；D 级是自己写进去的冗余；
   C 级是观测类别问题不是数据量问题）
8. **会犯什么错**（层间纪律 + 七条常见质疑回应 + 剩余缺口表）
9. **边界 + 一次翻车**（内部极小判据**结构性失明**；补规矩"每写一个判据必须写清它何时失效"）
10. 结尾（河与河床的比喻）+ 四级标记约定

**四级标记改为大白话**：【算出来的】/【我假设的】/【只在模型里算过】/【还没验】。

**不能声称的**（同 `16`）：A_D 未获现实证据；C 级读不出；数值实验只是排除程序写错。

---

## v1.6.9 — 2026-10-23

### XS 框架 v3.2 规范固定版（重构稿）：三层 + 一条定律

> **本轮不改内容，只改结构。**
> 用户指令：「基于现有的工程验证，我们重新来构架 XS 框架。」
> 审计 `14` 已把补丁归因到源头；本轮据此**重写整体骨架**，把 v3.1 的
> 「主文 + 九篇补篇」树状累积，改写为 **三层（L0 构成 / L1 约束 / L2 动力学）+ 一条定律**的柱状结构。

#### docs: `docs/xs-framework/16-v3.2-self-consistent-theory.md`（**可分享的理论文章**，现行主线）

**第二份产出，面向「读起来是一套理论」而非「读起来是一份台账」**：
不删内容、不加应用，只做三件事——**重排论证、分离理论与证据、划清边界**。

- **四级标记贯穿全文**：【定理】纯逻辑推出 ｜【猜想】未验证的关键假设 ｜【原型确认】仅原型内部算过 ｜【待验证】有程序未执行。
- **论证链**：L0 生成（三分 + 偏离→归宿矩阵）→ L1 约束（双约束、`P/Q/U`、混合偏倚）→ **§4 瓶颈：U 与 Ψ 之间没有桥** → §5 记法分层（S9：存量 vs 通量）→ §6 猜想 A_D 与它的两处问题（λ-冗余、单一驱动）→ §7 规范固定与主定理 → §8 五个投影 → §9 楔形分类 → §10 可辨识性四级 → §11 需要哪一类观测 → §12 猜想登记表 → §13 边界与自洽性三重检验。
- **§12 猜想登记表**：每条未闭合项写明「陈述 / 承担了哪些结论 / 证否程序 / **失败后果**」。
- **§13.3 自洽性三重检验**：逻辑闭合、参数闭合（`(E,U)` 给 2 自由度 ↔ A+B 级 = 2，与独立路径吻合）、反例闭合。
- 明确**不声称**：A_D 已证；B-B22 已闭合；共线性法在任何数据上优于峰距法（噪声方差对比未做）；R43′/R32/R45 已验证；本文是关于现实的理论。

#### figs: `figs/v32_wedge_classification.svg`（图 1）

`(s, ρ)` 平面楔形分类：竖线 `s=1` 与对角线 `ρ=s` 划出四区；
两片楔形以 `(1,1)` 为顶点，张开度 `|1−ρ|`；轨迹沿 `s` 下降方向穿过。

#### docs: `docs/xs-framework/15-v3.2-gauge-fixed-framework.md`（重构对照与工程台账）

**核心主张**：补篇二/五/六/七/八/九**六篇**（占理论篇幅约三分之二）**全部是同一条定律的五个投影**：

$$\frac{\mathrm dU}{\mathrm dt}=\rho_{exp}\cdot\frac{\mathrm dE}{\mathrm dt}+\kappa_E\,\rho_{exp}\,(1-\rho)$$

| 新增 | 内容 |
|---|---|
| **D1 层间纪律** | 归宿（收敛/发散/震荡）在 L0–L1，形状（峰序/峰距/双相）在 L2；**禁越层读取** ⟹ 一次性解释 B-B34 的 36.8%–76.5% 偏差带 |
| **U5 规范不变性判据** | 新增元规则：真值随 λ 改变的命题**不予登记**为可检验推论（预防 T17 的规范冗余型伪不可辨识） |
| **记号分层（沿用 B-B19，非新发现）** | `Ψ` 与 `U` 被写成同一式 `P − Q` 却分属通量层与存量层 ⟹ **解释了 B-B19 为何长期不闭合** ⟹ 桥即仿射定律。⚠ **作者更正（2026-10-23 晚）**：初稿称此为「新发现的第九个病根 S9」，实为 **v3.1 补篇四 §2.4 已登记**；v3.2 的贡献是①可执行分层规则 ②发现两处遗留旧记号未加取代标注 ③B-B19 获得闭形式与可证伪形式 |
| **命题 5.2 楔形分类** | `(s,ρ)` 平面被 `s=1` 与 `ρ=s` 两线划分：**钻研与反刍是同一楔的两半**；中间相 ⟺ `ρ ≠ 1`；张开度 = `|1−ρ|`。定理 5.2′ 由此从 4 参数压缩为 2 参数 |
| **可辨识性四级分类** | A `ρ_exp`（自足）/ B `sign(1−ρ)`（自足）/ C `ρ`·`κ_E`（需外部单位）/ D 四系数绝对值（**伪目标**）；A+B = 2 与补篇九定理 6.12 独立吻合 |
| **命题 6.1 第二可辨识视界** | 纯 `(E,U)` 观测下 `(κ_E, ρ)` **不可分离** ⟹ 审计的 R42 实为**循环** ⟹ 新缺口 **B-B45** |
| **命题 6.3 实验性规范固定** | λ-规范自由度源于 `a,b` 的**单位未固定** ⟹ **固定动作计数**即可测出 `(e₊,e₋,u₊,u₋)` ⟹ C 级闭合（**R45**） |
| **R43** | 共线残差 ≡ `g(U)` 沿路径的变异 ⟹ 残差是**信号**不是噪声（B-B43 部分推进） |
| **R43′** | 桥（B-B19）的可证伪形式 `β₀ = −ρ_exp·Ψ(k*_U)`，需 `Ψ` **独立标定**（用 `ΔE` 代替则退化为重言，本稿自纠） |

**两处修正审计 `14`**：

1. **R42 不是标定路径**（循环）⟹ 降级为记账恒等式；缺口改记 **B-B45**。
2. **B-B22「A_D 无任何有效证否途径」不准确** ⟹ S2 落地后陷阱 T16 不再构成假阳性 ⟹
   **补篇九 推论 6.13 本身就是 A_D 的有效证否途径** ⟹ B-B22 的真正内容收缩为
   「A_D 的**经验地位**」，严重程度由「框架自证困难」下调为「需外部经验证据」。

**缺口收敛**：真理论缺口 **21 项 → 4 项**（B-B22 / B-B42 / B-B43 / B-B44）+ 新增 B-B45。

#### docs: `docs/xs-framework/README.md`

新增 `15` 行（B-B22 重新定性、B-B19 闭合、B-B42 改写、B-B45 新增、读法建议 5–6 条）。

#### 不能声称的

A_D 未被证明；B-B22 未闭合；R43 / R43′ / R45 **均未做实验**；R30 结论仍限于四个调度族。

---

## v1.6.8 — 2026-10-22

### 工程补丁的源头归因审计（兼 v3.2 重构提案）+ R30 实验

> **本轮不是补一个洞，而是回头问：这些洞是不是我们自己挖的。**
> 用户判断："很多补丁在工程上很难完成，但源头上的一点修改往往可以一劳永逸。"
> 本轮据此把补篇一～九、实验 R18–R29 留下的**全部 21 项未闭合项**向上归因到源头的
> 定义/参数化，区分「真缺口」与「由源头选择不当造成的伪缺口」。

**结论：21 项 → 8 个源头病根 → 消失 7 / 降级 6 / 残留 7（其中真理论缺口 4）。**

#### docs: `docs/xs-framework/14-v3.2-source-refactor-audit.md`（审计）

八个源头病根（S1–S8），每个给出「症状 → 病根 → 修正 → 消灭哪些补丁」：

- **S1 参数化冗余 ★★★** — 补篇五 §1.2 的 `dE = a·e₊ − b·e₋` 把 **6 个量写成了 3 个自由度的模型**；
  多出来的就是定理 6.12 的 λ-规范自由度。规范固定后（`s:=r/θ_E`、`κ_E:=b·e₋`）内禀参数降为
  **3 个无量纲数 + 1 个时间单位**，λ 自由度**在定义层被消灭**。
- **S2 单一驱动是定义成分 ★★** — 定理 6.10/6.11 把"共享 r"列为**额外前提**，但 A_D 原文
  （"同一个 a、同一个 b"）已含它。写入**单一驱动附款**后，恒等式变成 A_D 的**代数推论**。
- **S3 换可观测量 ★★★** — 峰距法读的是 ρ 的**数值**（需要 κ_E ⟹ 调度依赖）；
  共线性读的是 ρ_exp 与 ρ 的**符号**（不需要 κ_E ⟹ 调度无关）。
- **S4 尺度算子内置** — 用补篇二已有的 `Ψ = ψ̄ + ξ` 定义分离尺度 τ\*。
- **S5 `m_eff := 2^U`** — U 是熵、m 是计数；桥由**定义**给出，不再是待证命题。
- **S6 状态依赖是必需的 ★★** — `U ≥ 0` 是定理而线性 A_D 会跑出状态空间 ⟹ **要求** A_D⁺ 的
  `g(U)`（`g(0)=0`）。⟹ 工程上被当作"违反族"的**状态依赖型，其实是 A_D 的正确形式**。
- **S7 偏离→归宿改矩阵** — 原 1:1 映射的**元数错了**，实际多对多。
- **S8 混合偏倚 = `I(M;Ω|φ)`** — 这是 `P := H(M|φ,Ω)` 定义的**直接推论**，不是未知量。

#### experiments: `experiments/r30_gauge_fixing/`（R30，仅标准库，三项验证）

| 验证 | 设计 | 结果 |
|---|---|---|
| **V1** 仿射恒等式逐样本成立 | 3×θ_E × 5×ρ × 4×调度族 = 60 组 | 最大残差 **1.998e-15**（相对 2.26e-16） |
| **V2** λ 规范不变性 | λ ∈ [0.25, 4]（16×） | θ_E 变 16×，而 `max\|dE\|=max\|dU\|=max\|ds\|= ` **0.0e+00**；(ρ_exp, ρ, κ_E) 恒定 |
| **V3** 调度无关性 | 同一批轨迹比较共线性法与峰距法 | 共线性 **0.000e+00**（全族恒定）；峰距 **0.38% / 63.02% / 57.73% / 20.75%** |

#### 核心恒等式（本轮最重要的单条结果）

$$\frac{\mathrm{d}U}{\mathrm{d}t} \;=\; \rho_{exp}\cdot\frac{\mathrm{d}E}{\mathrm{d}t} \;+\; \kappa_E\,\rho_{exp}\,(1-\rho)$$

**补篇六峰序、补篇八守恒、补篇九共线性三层技术装置 = 这一条关系的五个投影：**

1. 共线性 ⟹ 残差判据（定理 6.11）
2. 斜率 = ρ_exp ⟹ ρ_exp 可标定（**闭合 B-B23 一半**）
3. 截距符号 = sign(1−ρ) ⟹ 峰序（推论 6.4，只需符号，不需 κ_E）
4. 积分 ⟹ `U − ρ_exp·E` 仿射守恒（定理 6.10）
5. 截距**绝对值** ⟹ 需要 κ_E ⟹ **这正是峰距法调度依赖的根源（B-B34 被解释，非被拟合）**

#### 补丁处置

- **消失 7**：B-B23（主体）、B-B34、B-B35、B-B38、B-B40、B-B36、B-B11
- **降级 6**：B-B25（→需互信息估计器）、B-B37（→操作问题）、B-B26（→偏离多大）、
  B-B31（→η 原则）、B-B27′（→相对于 τ\*）、B-B28（→m_eff 连续）
- **缩减 1**：B-B39（状态依赖族变为必需项，剩阈值/时滞两族）
- **残留 7**：B-B22（A_D 经验地位，**唯一源头改不动的**）、B-B19、B-B29/32/33、B-B30、
  B-B24/B-B17/B-B12（纯工程）

#### 新缺口 / 新推论

- **B-B42（G1）**：κ_E 的外部标定 —— B-B23 的真正剩余。**已给出闭合路径**：
  由 U 峰处 `Ψ(k*_U) = κ_E(ρ−1)` 得 `κ_E = Ψ(k*_U)/(ρ−1)`（**R42**），走**补篇四已有的
  `Ψ = P − Q` 标定**，而非不可观测的 r。⟹ **把 B-B23 从"卡在不可观测的 r 上"改为
  "卡在可观测但需标定的 Ψ 上"——这是实质推进。**
- **R41（新推论，可检验）**：`β₀ = −ρ_exp·Ψ(k*_U)`，即**仿射截距 = −ρ_exp × U 峰处的势**。
- **B-B43（G3）**：A_D⁺ 的 `g(U)` 形式未知 ⟹ 补篇九全部结论须加"U 远离 0"的辖域限定。
- **B-B44（G2）**：尺度分离阈值 η 的取值原则未定。

#### 本轮两个新方法学陷阱

- **T17** ★：两种"不可辨识"必须分开 —— (a) **规范冗余**造成的**伪**不可辨识（改参数化即消灭）；
  (b) **真**不可辨识（如 κ_E 需要外部单位）。**把 (a) 当 (b) 会徒劳地去找外部标定。
  B-B23 的主体属于 (a) 而非 (b)。**
- **T18**：调度依赖性的来源**可定位到具体参数** —— ρ 的**符号**不需要 κ_E，ρ 的**数值**需要。

#### 不能声称的（审计 §8 摘录）

1. **不能说 A_D 被证明了** —— 规范固定只让它更干净、更可检，没给它经验地位。
2. **不能说 B-B23 完全闭合** —— 闭合的是 ρ_exp 与 ρ 的符号；ρ 的**数值**仍需 κ_E。
3. **不能说共线性法在任何数据上都优于峰距法** —— R30 是**无噪声**验证，有噪声时方差对比未做。
4. **R41 尚未验证** —— R30 未测 Ψ（只用 A_D 内部一致性），需要独立测 Ψ 的实验。

#### tests

新增 `tests/test_r30_gauge_fixing.py`（17 项，锁结构不锁数字）。
R 系列全套 **109 项通过**；仓库全套 **226 → 243**。

---

## v1.6.7 — 2026-10-20

### XS 补篇九（A_D 证否签名的判别力分解：守恒量、规范自由度与操作规程）+ R24 实验

本次**不改动任何既有模块的公开接口**；既有 199 项测试全部保持通过，新增 27 项 → **全套 226 项**。

> **本轮最重要的结论是否定性的**：A_D 唯一已登记的证否途径（R24）经检验**无效**。
> 详见下方「R24 的失效」。

#### 来源

补篇五 §6 的 R24 登记与 B-B22 缺口。R24 被登记为 A_D **唯一**的已知证否途径，
但登记时只有一句话，未经检验。本轮补上这个检验，并给出取代方案。

#### docs: `docs/xs-framework/09-v3.1-supp9-ad-refutation-signature.md`（补篇九）

- **定理 6.11（任务内线性失配守恒）★**：在 `A_D ∧ 共享驱动` 下，对**每一条样本路径**
  （不取期望、不取极限）：
  $U(k) - \rho_{exp}E(k) = \alpha + \beta k$，其中 $\beta = dt\,e_-(\rho_{exp}-\rho_{res})$。
  等价地 $\mathcal{R} := \min_v \lVert\Pi^\perp(vU - E)\rVert \equiv 0$，且 $v^\* = 1/\rho_{exp}$。
  ⟹ 增量云 $\{(dE_k, dU_k)\}$ **共线**，斜率 $\rho_{exp}$、截距 $\beta$。
- **推论 6.11′**：R24 的象限签名是**共线性的推论**，因而严格弱于共线性。
  由共线性可算：ρ>1 ⟹ `(−,+)` 恒空；ρ<1 ⟹ `(+,−)` 恒空；ρ=1 ⟹ 两者皆空。
- **定理 6.12（规范自由度 / 标定不可辨识性）★**：变换
  $G_\lambda: r_k\mapsto\lambda r_k,\ e_+\mapsto e_+/\lambda,\ u_+\mapsto u_+/\lambda$（`e₋, u₋` 不变）
  保持 $(E(k), U(k))$ **逐点不变**。⟹ 四个系数的**绝对值在任何只依赖 (E, U) 的
  可观测量下都不可辨识**；独立规范不变量恰为 **2 个**：$\rho_{exp}, \rho_{res}$（$\rho=\rho_{res}/\rho_{exp}$）。
  **这是 B-B23 不闭合的真正原因——不是数据不足，是结构性的。**
- **推论 6.8（临界违反强度）**：二次注入下 $\varepsilon_c = e_+/(4\theta_E)$；
  $\varepsilon\ge\varepsilon_c$ ⟹ Ψ 无实根、E 全程下降、峰退到 k=0。
- **推论 6.13（取代 R24）★**：A_D 的有效证否判据是**任务内线性失配残差** $\mathcal{R}$；
  跨任务并存只能证否合取 `(A_D ∧ ρ-任务不变性)`，**不能**证否 A_D 单独。

#### experiments: `experiments/r24_ad_refutation/`（R24，仅标准库，九个臂）

- **臂 A（预言复核）**：15 组配置，`cross_task_consistent` **15/15**；
  `max|d_err| = 0.00148`（≤ 离散化界 0.00417）；`max resid = 1.22e-14`；
  `max nmin = 0`；`v*` 误差 4.0e-14 ⟹ P1–P4 对一般 `u₊` 成立。
- **臂 B（混淆臂，★核心）**：
  - η=0 ⟹ 跨任务判据误判率 **0.0000**（无误报）✓
  - η>0 ⟹ 误判率升至 **1.0000**（**假阳**）⟹ **T12 确认**
  - 同一批数据上，任务内残差判据恒 ~1.4e-14、任务内签名唯一性恒 0.0000 ⟹ **判据 2/3 正确**
- **臂 C / E（分离力与视界，负面结果 N3）**：内部极小法与任务内双象限法
  在**全部网格点**上探测率恒为 **0.00** ⟹ 结构性不可检（见下「为什么」）。
- **臂 E2（阈值标度，★负面结果 N4）**：翻转点**不是** `w=1`，
  且随 θ_E 变化（θ_E=0.5/1.0/2.0 ⟹ w\* ∈ (0.20,0.40] / (0.40,0.60] / (0.90,1.00]）
  ⟹「峰退到左边界」是**模型依赖**判据，不能独立使用。
- **臂 D（共享 vs 独立噪声）**：共享 ⟹ 残差恒 ~7e-15，**与 s 无关**；
  独立 ⟹ 残差 >0 且 ∝ s（0.0094 / 0.019 / 0.039 / 0.078）⟹ **T16 的来源**。
- **臂 F（SOP 端到端）★**：三分类 **15/15 全对**；同一批数据上 R24 原判据 **5/5 误判 D3**。
- **臂 G（四系数辨识）**：r 可观测 ⟹ 相对误差 ~2.3e-14（**精确**）；
  r 不可观测 ⟹ OLS 仍一致但抖动随 s 增长（最大 4.5e-2）。
- **臂 G2（规范自由度）★**：λ 跨 16 倍，四个系数全变，而 `(E(k), U(k))` 差值 **0.00e+00**。
- **臂 H（噪声地板与平均律）★**：ε=0 ⟹ log-log 斜率 **−0.511 / −0.490 ≈ −1/2**；
  ε=0.05 ⟹ 斜率 ≈ 0，收敛到系统地板 **0.0407 = 0.81459 × 0.05**。

#### R24 的失效

> R24（补篇五 §6）**证否的是合取 `(A_D ∧ ρ-任务不变性)`，不是 A_D 单独。**
>
> 「ρ 依赖任务」——A_D 在每个任务内部**严格成立**——会产生**完全同样的**并存现象。
>
> ⟹ **B-B22 保持开放且升级：A_D 目前无任何已登记的有效证否途径。**
>
> 取代方案：**任务内线性失配残差**（定理 6.11 / 推论 6.13），
> 它在 `A_D ∧ 共享驱动` 下**逐样本恒为 0**，且灵敏度系数与 θ_E 无关。

#### 本轮五个方法学陷阱 T12–T16

- **T12** ★：跨任务并存 ≠ 违反 A_D（ρ 任务依赖产生同样现象，实测误判率 → 1.0）。
- **T13**："是否报证否"是逐 rep 的随机变量，必须报率。
  （本轮实际踩到：臂 B 的 B3 判据因符号格式不一致而恒为 0，属**假阴性**，已修正。）
- **T14**：判别"内部极小"必须用显著度过滤，否则噪声下的假变号淹没信号（承 T5/T11）。
- **T15**：ε 过大时 Ψ 无实根 ⟹ 探测率关于 ε **非单调**（已量化，见臂 C/E）。
- **T16** ★：观测层独立噪声也会使残差非零 ⟹ 残差判据只证否 `(A_D ∧ 共享驱动)`。
  修正：**重复测量后平均**再算残差（臂 H 验证 `1/√N` 律）。

#### docs: 补篇五 §6 加 R24 取代标注 / B-B22 升级标注 / B-B23 定性标注（U1）；主文加标注 ⑩；缺口表续 B-B38–B-B40

---

## v1.6.6 — 2026-10-20

### XS 补篇八（峰距标定程序：可辨识性、视界与操作规程）+ R29 实验

本次**不改动任何既有模块的公开接口**；既有 184 项测试全部保持通过，新增 15 项 → **全套 199 项**。

#### 来源

补篇六 §3.6 推论 6.5 与 §5 的 R29 登记。补篇六给出了一条**公式**，
补篇八检验把它当作**估计程序**用时的边界——公式正确不等于程序可用。

#### docs: `docs/xs-framework/09-v3.1-supp8-peak-distance-calibration.md`（补篇八）

- **定理 6.10（峰序守恒，路径级）★**：在 A_D 下，若 E 与 U 由**同一** r 过程驱动，
  取 $U^s := U/\rho_{exp}$（不改变峰位），则
  $U^s(k)-E(k) = \text{const} + dt\,e_-\!\left(1-\frac{\rho_{res}}{\rho_{exp}}\right)\!(k+1)$，
  **与 r 及噪声完全无关**。⟹ 峰序**逐样本守恒**，判反率为 **0**
  （实测：所有 s 档、所有 ρ 档 err_rate 均为 0.0000）。
- **推论 6.5′（可辨识性的准确表述）**：R29 可辨识的是峰距 $\Delta$ 本身，
  **不是** $\lvert\ln\rho\rvert$ 的绝对值。后者还须**独立知道 r 调度族**。
  **这是 B-B23 只能"部分闭合"的根因。**
- **视界律**：$\ln\rho$ 的可测窗口**宽度恰为 $L=\lvert\ln(r_{hi}/r_{lo})\rvert$**，
  中心为 $-\ln\theta_E$。超界时程序**不报错**，只是测不到峰。
- **操作规程**：两条件步数 K 不同时**必须归一化** $\delta=\Delta/(K-1)$。

#### experiments: `experiments/r29_peak_distance_calibration/`（R29，仅标准库，六臂 + 恒等式臂）

- **臂 A（公式精确性）**：无噪声下峰距命中解析预测，最大残差 **0.898 步** ≈ 1 步的
  **离散化硬界**（峰位 = 最后一个增量为正的索引，与连续穿越点差至多 1 步）。
  灵敏度 $C=(K-1)/L=54.614$ 步 / e 倍；过原点拟合 $C_{fitted}=54.506$（偏 −0.20%）。
- **臂 B（恢复精度）**：估计量 $\widehat{\lvert\ln\rho\rvert}=\bar\Delta/C$ **无偏**
  （最大 bias ≤ 0.018，随 $n_{rep}$ 缩小 → 来自离散化，非估计量本身）；
  $sd(\widehat{\ln\rho})\approx 0.35\,s$。给出可操作的样本量判据。
- **臂 C（比率设计，C 未知）**：
  - **C1**（同 K=241）：相对误差 **0.0015–0.0273** ✓ 比值消去 C 成立。
  - **C2**（异 K：241 vs 121）：raw 比值相对误差 **0.79–1.01**（约 2 倍系统偏）；
    归一化后 **≤ 0.1057** ✓。⟹ **归一化消除系统因子，但不恢复分辨率损失**
    （K 减半 ⟹ Δ 减半 ⟹ 同样的 ±1 步离散误差相对翻倍）。
- **臂 D（调度误设，负面结果）**：真调度非对数斜坡时，按推论 6.5 反解的偏：
  log_ramp 0.039（基线）/ piecewise 0.368 / sigmoid 0.548 / **linear 0.765**。
- **臂 E（可辨识性视界）**：窗口 $\ln\rho\in(-2.1972,+2.1972)$，宽度 = L ✓；
  $\ln\rho=2.1$ 仍可测，$\ge2.2$ 峰落窗外、不可测 ✓。
- **臂 F（峰序分辨率）**：shared 下 err=0、tie 随 s↑ 而↑（退化为"不确定"）；
  indep 下 err 最大 **0.4125**（$s=0.8$，$\lvert\ln\rho\rvert=0.05$），随 $\lvert\ln\rho\rvert$ 增大快速降到 0。
- **臂 G（恒等式验证）**：$U^s-E$ 线性恒等式逐点残差 **5.684e-14**（全窗口 6.395e-14），
  含 $\rho_{exp}\in\{0.5,1,2\}$ 的一般情形 ⟹ **恒等式精确成立**。

#### 本轮四个方法学陷阱 T8–T11

- **T8**：不可跨 K 直接比峰距（$C$ 含 $K-1$，实测差因子 2.000）⟹ 必须归一化。
- **T9**：$\lvert\ln\rho\rvert$ 的**绝对值不可辨识**，只有比值可辨识（越界声称风险）。
- **T10**：峰靠边界时 argmax **向内偏**；须先判"峰是否为内点"，否则把不可测报成"ρ 比值很小"。
- **T11**：平滑窗口会**合并近邻双峰**，造成峰序判定假阴性 ⟹ 臂 F 用 $w=1$ + 显著度。

#### docs: 补篇六 §5 加 R29 结果指针（U1）；主文加标注 ⑨；补篇四缺口表续 B-B34–B-B37

#### fix: LaTeX 转义损伤修复（命中 5 行，零误伤）

对全量 md（工作区顶层 + 仓库 `docs/`，共 48 个文件）做「连续两个反斜杠」扫描：
命中 5 行，全部归一为单个反斜杠。脚本会先向前扫描最近的 `begin`/`end` 环境，
若该行处于 `aligned` / `cases` / `matrix` / `array` 等**合法换行环境**内则跳过；
本次 5 行均不在其中，故**无误伤**。

- **本轮新引入 3 处**：`docs/xs-framework/README.md` 第 28 行（定理 6.10 公式，多处 `\text` `\,` `\!` `\left` `\frac` 被写成双反斜杠）、第 89 行（B-B36 的 `$U\ge0$`）；顶层补篇四第 237 行（同）。
- **历史遗留 1 处（两份副本同步修复）**：`XS-框架论文-统一重构版.md` 与 `docs/xs-framework/10-v3.0-unified-framework-paper.md` 第 85 行，「极小可审计」一栏的 $B\setminus\{b\}$ 被误写成「`B` + 双反斜杠 + `{b}`」的形式。

修复后全库重扫：双反斜杠 **0** 处、非 UTF-8 **0** 个文件、控制字符 **0** 个。

**教训（记入纪律）**：本轮 3 处新损伤说明「经命令行写入含反斜杠的 LaTeX」仍是高风险路径；凡含 `$...$` 的长段落，一律**用文件写入工具或脚本文本块**落盘，并在每轮收尾**重扫一次双反斜杠**。

## v1.6.5 — 2026-10-19

### XS 补篇七（微观回落的本体性、收敛的局部性、峰的存在性）+ R28 实验 + 补篇六的再修正

本次**不改动任何既有模块的公开接口**；既有 169 项测试全部保持通过，新增 15 项 → **全套 184 项**。

#### 来源主张（用户四条）

1. "增 / 减"若基于统计，则**都可能存在微观回落**。
2. 有效认知 = 收敛，但**不保证排除伪收敛**。
3. 收敛可以是**过程的、局部的**，而全局可能不收敛。
4. 峰**理论上一定存在**（X：无 → 值 → 认知）；强调"有效认知"是因为它**保证收敛 ⟹ 保证减相**，
   而发散不一定有减的过程。

#### ⚠ 对补篇六 §1.3 的**再修正**（I 型，上一轮判断错误）

- **CE-2（纯净消解）撤回**。在"$m(0)=0$（从无）"的前提下，$m(0)=6$ 只能是完整过程的
  **后一段**；作为一段，它前面必有 $0\to6$ 的扩张，整条曲线有峰。
- **CE-1 保留但降级**为平凡过程（一步命中 $0\to1$），与命题 6.1 的"非平凡"同一排除项。
- **`m(0)\geq2$ 这个非平凡条件写错了**（在 $m(0)=0$ 下自动失效）。正确条件是
  **A_expand: $\exists k:\ m(k)\geq2$**，登记为 **B-B30**（是假设，非定理）。
- **"有效"的作用域被上一轮读错**：用户强调"有效"是因为**有效 ⟹ 收敛 ⟹ 有减相**
  （发散不保证），不是因为"有效 ⟹ 有峰"。
- **净结果：在"从无 + 非平凡"下峰必然存在**（命题 6.1′），比补篇六的结论**更强**——
  上一轮收缩过度了。

#### docs: `docs/xs-framework/08-v3.1-supp7-micro-reversal-locality.md`（补篇七）

- **命题 6.6（微观回落率）★**：A_D^stat 下 $P(dU<0)=\Phi((\ln\theta_U-\ln\bar r+s^2/2)/s)$，
  **精确式非近似**（因为 $dU<0\iff r<\theta_U$ 本身就是阈值事件）。
- **推论**：峰处回落率 $=\Phi(s/2)>1/2$，随 $s$ 单调增——直觉（"峰顶增量为 0 → 1/2"）是错的，
  因为对数噪声的中位数 $<$ 均值。
- **命题 6.7（尺度分离）**：单峰在**路径层是假命题**，在 **E 层是真命题**。
  **B-B27 升级为 B-B27′**（从"测量问题"更正为"过程的本体属性"）。
- **推论 6.7.1（操作规程）**：定义"峰"**优先用显著度 $h$，不要用平滑窗口 $w$**
  ——$h=0.1$ 在 $w=1$ 时就恢复单峰率 1.000，而平滑要 $w=41$ 且**会移动峰位**。
- **命题 6.8**：局部收敛 $\not\Rightarrow$ 全局收敛（$F$ 可继续添加肯定事实使新宏可加载）。
- **定理 6.9（峰数–回跳拓扑律）**：$\lvert n_+-n_-\rvert\leq1$；止于下降段 $\Rightarrow n_+=n_-+1$。
  **显著度过滤会破坏该恒等式**（两可比峰 + 一浅谷 → 谷被滤掉、两峰保留）。
- **重开 $\subseteq$ 回跳**：重开要求回落**确实到达过 $m=1$**。
  **无回跳 = 形状单峰**；**无重开 = 一次性收敛**；后者蕴含前者，反之不然。
- **推论 6.11**：震荡 = 局部收敛被反复重开（给了震荡一个可数的形状判据）。
- **推论 6.12**：收敛 = 正确收敛 + 伪收敛；R21 数据显示伪收敛占比随内约束
  **12.4% → 37.6% → 53.1%**：**剩下的"有效认知"越来越不可信**。

#### experiments: `experiments/r28_micro_reversal/`（R28，仅标准库）

- **关键建模改动**：把涨落从**观测层**移到**过程层**（$r=\bar r\,e^{sZ-s^2/2}$）。
- A1 解析律命中：$s\in\{0.1,0.2,0.4,0.8\}$ 下全轨迹 max|偏差| $\leq0.0225$，
  MC 标准误 0.0079，**全部在 3se 内**；峰处回落率 0.514/0.529/0.577/0.658 vs $\Phi(s/2)$ 0.520/0.540/0.579/0.655。
- A2 **路径层单峰率 0.000**（$s\geq0.2$）→ $w=41$ 恢复到 0.843–1.000；峰数在 $w$、$h$ 上均单调不增。
- 臂 B：三段非单调调度，峰 $=[19,179]$，回跳 1，**恒等式 2=1+1 成立**；噪声下仍为 2.000。
- 臂 C **不依赖 A_D** 的生灭模型：扩张相**内生**微观回落率 0.134（模型无任何外加噪声）；
  **局部收敛 100% vs 全局收敛 13–34%**（扫 $\tau_c\in\{1,3,5,10\}$）；
  交替律在 **100%** 路径上成立。

#### 本轮三个方法学陷阱 T5–T7

- **T5（假阳性，性质最严重）**：用按**观测噪声**标定的容差判单峰，而噪声在**增量**上 →
  随机游走累积偏差远超每步容差，**单峰率假性 100%**。修正：主判据改为不依赖容差的"显著峰数 == 1"。
- **T6**：假设"峰数 = 回跳数 + 1"在显著度过滤后仍成立 → 实测仅 0.458。
- **T7**：用"止于上升/下降"解释残差 → 条件化后反而更低（0.301）。真因是 T6 的过滤。

#### docs: 补篇六加三处 U1 交叉引用标注（§1.3 CE-2 撤回、§2.3 回跳精确化、§3.2 非平凡条件、§4.3 B-B27 定位更正）

## v1.6.4 — 2026-10-19

### XS 补篇六（认知过程的双相形状：X 单调 / M 单峰）+ R27 实验 + Markdown 转义损坏修复

本次**不改动任何既有模块的公开接口**；既有 155 项测试全部保持通过，新增 14 项 → **全套 169 项**。

#### 来源主张

用户提出：完整有效的认知过程不仅在前期认知事实（X）的递增，也在认知（M）的迭代——
**增到峰值时，又有一个收敛的减的过程**。

#### docs: `docs/xs-framework/07-v3.1-supp6-biphasic-shape.md`（补篇六）

- **承载量选择**："M 的增/减"操作化为 $m(k) := \lvert A_k\rvert$（与主文三分定义逐字对齐：
  收敛 ⟺ ∃! M\*），U 为辅。
- **一处必须修正（E 型）**：原话的"完整**有效**"须拆开。**"有峰"不是有效性的条件**——
  CE-1（直接命中）与 CE-2（纯净消解）都是有效但无峰的反例。
  **"有峰"是完整性的标志。** 修正后主张才可证伪。
- **与既有结果的关系（诚实）**：补篇一极性定理 3.3 已蕴含"混合极性 → 涨落"，
  补篇一 §3.8 的三层分工已说明"减"被移位到 eff 与 $A_k$。本补篇新在两处：
  **N1 相序**（要求"肯定相在前、否定相在后"，强于"混合"）与
  **N2 峰数 = 回跳数 + 1**（单峰 ⟺ 无回跳）。
- **命题 6.1**（**无条件**）：收敛且非平凡 ⟹ 存在消解相。只用主文定义 5.1，**不依赖 A_D**。
- **定理 6.3（峰位）**：$k^\*$ 由 $r(k^\*)=\theta_U$ 给出，即**探索-利用切换点**，与噪声无关。
- **推论 6.4（峰序定理）★**：$\rho_{res} > \rho_{exp} \iff k^\*_U < k^\*_E \iff$ 中期为钻研相 $(+,-)$；
  反之为空洞相 $(-,+)$。**把定理 5.2′ 推进到时间域**。
- **推论 6.5（峰距标定式）**：$\lvert k^\*_E - k^\*_U\rvert \propto \lvert\ln(\rho_{res}/\rho_{exp})\rvert$
  —— **B-B23 的第一条标定途径（部分闭合）**。

#### experiments: `experiments/r27_biphasic/`（R27，仅标准库）

- 对数斜坡 $r: 9 \to 1/9$，$K=241$，$n=400$/档；两配置（$\rho_{res}=3$ 与 $\rho_{res}=1/9$）。
- **峰位命中解析预测**（平台比对，避免并列取小索引的伪像）；
  **峰序吻合率 1.000**；中期相签名与推论 6.4 逐项一致；对照臂确认"残缺过程无峰"。
- **负面结果一（B-B27，高）**：单峰性对观测噪声**不鲁棒**——未平滑时检出率仅 0.045
  （$\sigma=0.15$），$w\ge5$ 后恢复至 $\ge0.985$。**必须在声明的粗粒化尺度上检验。**
- **负面结果二（B-B28）**：峰位估计依赖峰后陡度；B 组的 184.5 只能读作"峰位 $\ge180$"。
- **踩到并修正的四个方法学陷阱**（均写进脚本 docstring 与报告）：
  T1 用内置 `hash()` 做种子（str hash 按进程随机化 → 改用 `zlib.crc32`）；
  T2 峰位点值比对（并列取小索引的伪像 → 改用平台比对）；
  T3 用平台下界估峰位（噪声下**系统性下偏** 42 vs 60 → 改用平滑序列 argmax）；
  T4 平滑后仍用固定容差（→ 容差取 $3\sigma/\sqrt w$）。

#### fix: 5 个 Markdown 文件的 LaTeX 转义损坏（36 处）

- **发现**：`docs/xs-framework/README.md`、`01-`、`04-`，以及顶层
  `XS-认知框架-二分-宏加载-概率三分.md`、`XS-认知框架-v3.1补篇三-…md` 中，
  LaTeX 宏被当成转义序列吃掉反斜杠：`\rho` → CR+`ho`、`\varphi` → VT+`arphi`、
  `\approx` → BEL+`pprox`、`\text`/`\tau` → TAB+`ext`/`au`。
- **修复**：逆映射精确可逆（控制字节 → 反斜杠 + 控制名首字母），共 36 处；
  已确认所有 TAB 均在 LaTeX 内部、无一是缩进。修复后全文扫描**残留 0**。

## v1.6.3 — 2026-10-19

### XS 补篇五（U/Psi 桥、内约束强度、P 的理解者依赖）+ R19 / R21 实验 + Lean 候选命题

本次**不改动任何既有模块的公开接口**；既有 141 项测试全部保持通过，新增 14 项 → **全套 155 项**。

#### docs: `docs/xs-framework/06-v3.1-supp5-bridge-constraint-strength.md`（补篇五）

处理四项缺口，按优先级：**B-B19 → B-B21 → B-B17 → B-B16**。

**① B-B19（高 → 中）：U 与 Psi 的桥 —— 已定位，未闭合**

- **猜想 5.2 作为无条件命题被证否**（登记 X−）。三个反例模型并存：
  正常收敛 $(-,-)$、越想越乱 $(+,+)$、钻研型 $(+,-)$。
  此前写"待证"，实际上**待的是否证**。
- 引入显式假设 **A_D**（展开 / 消解双速率线性分解）：
  $\mathrm{d}E = a e_+ - b e_-$，$\mathrm{d}U = a u_+ - b u_-$。
- **定理 5.2′**（取代猜想 5.2）：$\operatorname{sign}\bar\psi = -\operatorname{sign}(\mathrm{d}U/\mathrm{d}k)$
  $\iff \rho_{\text{res}} > \rho_{\text{exp}} \wedge r \in (\theta_E,\theta_U)$，
  其中 $\rho_{\text{exp}} = u_+/e_+$（展开的效率比）、$\rho_{\text{res}} = u_-/e_-$（消解的效率比）。
  **符号方向没错，错的是"无条件"。**
- **跨层引用规则**由补篇四的无条件禁令改为**条件化许可**：
  仅在已验证 $\rho_{\text{res}} > \rho_{\text{exp}}$ 时可互引（验证方法见 R23）。
- **推论 5.3（新结果）**：震荡（$\bar\psi = 0$）期内 $U$ 的走向由 $\rho_{\text{res}}/\rho_{\text{exp}}$ 决定
  —— $U$ 降 = **酝酿**，$U$ 升 = **反刍**，$U$ 不变 = 纯空转。
  **这使"有价值的酝酿"与"无价值的反刍"第一次可判别，且不需要等终局。**

**② B-B21（已闭合）：内约束强度的操作化**

- 撤回的 R18-b 有了替代：内约束强度 = **$O_{\text{adm}}$ 剪枝量**
  $Q = \log_2(|O_{\text{all}}|/|O_{\text{adm}}|)$。不依赖 $\tau$ 或注入策略，无同义反复。
- **同时证否了我自己提出的一条预言**（登记 X−）：
  "虚假唯一性"——预期剪枝会通过人为排除第二名来抬高伪收敛；
  实测条件于真值未被剪时 F2 **反而单调下降**（9.9% → 0.0%）。
- 该证否**修正了主文 §5.4 伪收敛成因 2**（E 型）：
  **过早锁定的危害条件不是"候选太少"，而是"正确选项不在 $A_k$ 中"。**
  候选少是**放大器**，不是**成因**。

**④ B-B16（已闭合）**：$P := H(\mathcal{M}\mid\varphi,\ \Omega)$（$\Omega$ = 理解者）。
由此：**三态不是句子的属性，是"句子 × 理解者"配对的属性**；
连带**新增 B-B25**（跨理解者混合估计系统性低估 $P$，偏倚 $= I(\mathcal{M};\Omega\mid\varphi)$，未给修正程序）。

#### feat: `experiments/r19_bridge/`（R19）

A_D 的**内部一致性检验**——不是经验验证（无人造参数以外的经验来源；
上一轮 R18 出现过"引用了不存在的数据"，本轮对策是凡写入结论的数字都来自本目录脚本的直接输出）。

| 配置 | $\rho_{\text{exp}}$ | $\rho_{\text{res}}$ | 预测可达 | 实测 | 一致 |
| --- | --- | --- | --- | --- | --- |
| I | 0.30 | 0.90 | $(+,+)\ (+,-)\ (-,-)$ | 同 | ✓ |
| II | 0.90 | 0.30 | $(+,+)\ (-,+)\ (-,-)$ | 同 | ✓ |
| III | 0.50 | 0.50 | $(+,+)\ (-,-)$ | 同 | ✓ |
| IV（违反 A_D：加 $-0.2a^2$） | 0.30 | 0.90 | 四种全部 | 四种全部 | ✓ |

- **配置 IV 是最重要的一项**：它证明定理 5.2′ 的**证否条件不是空的**——
  违反 A_D 线性性后，$\Psi = 0$ 从一个点变两个点，两种反号组合**并存**。
  若"并存"在原理上不可能，R24 这句证否就是废话。
- 噪声鲁棒性：$\xi \in \{0, 0.05, 0.15, 0.40, 0.80\}$ 五档，可达组合集合在每档都与解析预言一致。
- 推论 5.3 探针：配置 I 震荡期 $U$ 触底 100%（在震荡中完成确定 = 酝酿）；
  配置 II 触顶 100%（反刍）；配置 III 不变。
- **诚实处置**：非线性配置的"符合推论 5.3"字段返回 `None`（不适用），**不返回 True**，
  由 `test_nonlinear_is_marked_not_applicable` 锁死，防止假的确认。

#### feat: `experiments/r21_constraint_specificity/`（R21）

| $K=\lvert O_{\text{adm}}\rvert$ | $Q$ | 真收敛% | **F2%** | F1% | gap | $M^*$被剪% | 存活 F2% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.000 | 70.0 | 9.9 | 20.1 | 0.257 | 0.0 | 9.9 |
| 14 | 0.515 | 51.6 | 20.4 | 28.0 | 0.297 | 30.6 | 7.8 |
| **10** | **1.000** | 39.7 | **23.9** | 36.5 | 0.426 | 50.4 | 5.8 |
| 6 | 1.737 | 23.9 | 19.3 | 56.8 | 0.629 | 70.5 | 2.5 |
| 2 | 3.322 | 9.9 | 11.2 | 78.9 | 0.589 | 88.4 | 0.0 |

- **F2 倒 U**，峰值 $Q \approx 1.0$。10 个数据点，非单调性明确
  （4 档时不足以支撑该主张，故另写 `r21_sweep.py` 做细粒度扫描）。
- **两个效应必须分层**，否则会把"真值被剪"的**平凡效应**写成发现：
  效应 1 排错（正面，存活子样本真收敛 70.0 → 85.1、F2 9.9 → 0.0）；
  效应 2 排对（负面，$P(M^*\text{被剪}) = 1-K/20$，设计为必然）；
  效应 3 清空（负面，$A_k = \varnothing \Rightarrow$ 发散）。
- **高分辨率复核**（`r21_highres.py`，$N_{\text{FEAT}}=16$ / $N_{\text{MACRO}}=40$）：
  倒 U **保持**，峰值仍在 $Q \approx 1.0$ —— 排除"$\sigma$ 量化伪像"。
- **负面结果（登记 B-B24）**：按语义簇剪枝与随机剪同等数量**几乎完全重合**
  （cluster 70.4/54.9/37.8/18.3 vs random 71.2/54.4/38.4/18.5）。
  → 测到的全部效应都是**数量效应**；原型缺乏语义结构，
  → **该操作化应称为「$O_{\text{adm}}$ 剪枝量」，不应称为「措辞特异性」**。

#### feat: `formal/candidates/DualConstraint.lean`（B-B17，⚠️ 未编译验证）

- R18 的命题规格：DC1–DC8（结构性公理）+ T1–T6（可导出定理）+ 不可达登记 X1/X2。
- **放在 `candidates/` 是故意的**：该目录不在 `lakefile.toml` 的 `Formal` 目标下，
  **不会被 `lake build` 编译**。本机无 Lean 工具链（`lean` / `lake` 均不可用），
  未编译验证的文件进 `Formal/` 会直接破坏构建。
- 已做**静态检查**（声明顺序 / 前向引用 / 括号平衡 / namespace 匹配），
  修掉一处前向引用（`LoadFnPreserved` 在 DC1 之后声明）。但**静态检查不等于编译通过**。
- **形式化范围受限**：只编码不依赖 B-B19 的部分（两类失效的可分离性与签名 R18-a）。
  R18 的**定量签名**（$|\Delta\sigma| \le 0.006$ 且内外分离度相差一个数量级）
  是数值不等式，Horn skeleton **结构性不可表达**，登记为 X1，交仿真层回归测试锁死。
  「收敛 ⟺ $U \approx 0$」「定理 5.2′」「推论 5.3」因 B-B19 未闭合暂不编码，登记为 X2。
- `formal/GAPS.md` 新增条目 **C1**，含接入步骤。

#### test: 14 项新增（`test_r19_bridge.py` 7 项 + `test_r21_constraint_specificity.py` 7 项）

沿用 R18 的设计原则：**只锁结构性预言，不锁具体数值。**

- R19：两种反号组合的**互斥性**（由 $\rho_{\text{res}}$ vs $\rho_{\text{exp}}$ 决定）；
  退化配置两种都不可达；**非线性配置必须两种并存**（锁证否条件的判别力）；
  推论 5.3 的 $U$ 方向；非线性配置不得被判为"符合"。
- R21：F2 倒 U 且峰值在内部；存活子样本 F2 下降 / 真收敛上升；
  全样本真收敛下降 / F1 上升；$M^*$被剪率符合设计值；
  **倒 U 在高分辨率下保持**。
- **防回退（R21 最关键的一条）**：`test_source_does_not_touch_stop_rule_or_thresholds`
  检查源码**没有**缩放阈值、且判定函数唯一。
  R18-b 的失败正是用接受阈值编码"内约束强度"（同义反复），此测试防止重蹈。

#### 理论侧缺口（详见 `docs/xs-framework/README.md` 第四节）

未闭合（高）：**B-B22**（A_D 本身是未证的建模假设；证否签名 R24 未跑）、B-B11。
未闭合（中）：**B-B23**（$\rho$ 的单独标定未给）、**B-B24**（原型无语义结构）、
**B-B25**（$\hat P$ 混合估计偏倚未修正）、B-B16 已闭合、B-B17 部分推进。
已闭合：B-B6 / B-B10（主要）/ B-B13 / B-B14 / B-B21 / B-B16。
**B-B19 由「未证」变为「已定位未闭合」，优先级 高 → 中。**

> **注意**：B-B19 的处置**不是**闭合。断裂没有消失，只是被定位了——
> 它从"两个不可通约的层"变成"一个可测量的不等式"。
> 代价是引入 A_D（B-B22），收益是 R23/R24 与推论 5.3。

## v1.6.2 — 2026-10-19

### XS 框架理论文档集入库 + R18 双约束失效实验

本次**不改动任何既有模块的公开接口**；既有 134 项测试全部保持通过，新增 7 项 → **全套 141 项**。
版本号的说明见 README「已知局限」第 8 条：**XS 理论版本号（v3.0 / v3.1）与仓库版本号是两条独立的线**，
文档入库不触发仓库主版本升级；此处升 v1.6.2 是因为新增了实验代码与测试。

#### docs: XS 框架理论文档集（15 篇）

`docs/xs-framework/`，索引见该目录 `README.md`。按 **U1（取代不删除）** 纪律维护：
后出的补篇不修改前篇正文，只在前篇顶部插入条件化标注互指，故"同一命题有两个版本"是预期的，
**以编号最大者为准**。

- **v3.1 认知框架主线（5 篇，本轮新增）**
  - `01-…main-dichotomy-macro-loading-trichotomy.md` 主文：认知二分（X = 象 / S = 运算）
    → 事实公理 F 与运算公理 O → 宏加载（加载度 σ）→ 三分（收敛 / 发散 / 震荡）→ 概率三分。
    核心结果：**三条偏离 → 三个归宿**（事实集开放 → 概率性；资源有界 → 发散；承诺可撤销 → 伪收敛）。
  - `02-…supp1-chaos-boundary-negation-closure.md`：混沌边界与否定闭合。定理 3.2
    （**非单调性被局部化到 `eff` 这一个映射上**）；极性定理（**混合极性 = 震荡的动力学来源**）。
  - `03-…supp2-potential-drift-diffusion.md`：**势 Ψ := P − Q**，分解为 **Ψ = ψ̄ + ξ**
    （统计定向 + 局部回落）；定理 6.2″（ψ̄ ≠ 0 ⇒ 几乎必然触及边界，混沌被排除）。
  - `04-…supp3-dual-constraint-PQ.md`：双约束构成（外 = 输入物理约束 / 构成性 / 对 X；
    内 = 必须"理解" / 规范性 / 对 S）+ **两类失效分工**（外失效 → F2，内失效 → F1）。
  - `05-…supp4-calibration-commensurability.md`：**Q := I(M; C_in | φ)**，**U := P − Q = H(M | φ, C_in)**；
    P 的标定程序（经验分布熵 + Miller–Madow 修正）；三组二分**只按弱形式保留，不主张严格同构**。
- **v3.0 统一版与元规则（4 篇）**：统一重构版论文、三版统一与工程反哺、U1–U4 更新纪律 + 推论 P、流形裁决。
- **8 月原始稿与审计（6 篇）**：六篇去重正文、三项审计（向量空间隐喻 / 线性-非线性尺度进化 / 非专业观察者）。

#### feat: R18 双约束失效可分离性实验

`experiments/r18_dual_constraint/`（5 脚本 + 5 个 json + `report.md`），**仅依赖标准库，不需要 numpy**。
四臂只改各自那一个约束，不动落定规则 / 阈值 / 判定函数，F1 与 F2 由 σ 的演化路径**导出**而非指定。

- findings（主实验，800 次/臂）：

  | 臂 | 真收敛 % | F2 % | F1 % | 终局 σ |
  | --- | --- | --- | --- | --- |
  | base | 71.1 | 9.6 | 19.2 | 0.646 |
  | A 外失效 | 0.5 | **79.8** | 19.8 | **0.649** |
  | B 内失效 | 0.0 | 0.0 | **100.0** | **0.375** |

  **判别签名 R18-a 成立**：`|σ(A) − σ(base)| = 0.003`，而 `σ(base) − σ(B) = 0.271`。
  → 外约束是**方向约束**（错了但一样确信），内约束是**强度约束**（不够确信就停）。
- findings（剂量-反应，600 次/档）：外约束失效程度 `p_swap: 0 → 1` 时
  F2 率 **10.8 → 80.2 单调上升**，而终局 σ **保持平坦**（0.608–0.648，无趋势）。
  → 外约束失效的标志不是"想得不认真"，而是"**认真地想错了**"。
- findings（鲁棒性）：6 组 `(N_MACRO, RHO_OBS, P_NOISE, B)` 下签名全部成立，
  `|σ(A) − σ(base)| ≤ 0.006`，`σ(base) − σ(B) ≥ 0.199`，**分离度相差一个数量级以上**。
- findings（高分辨率复核 `r18_highres.py`，`N_FEAT=16 / N_MACRO=40`）：
  σ(base) = 0.587、σ(A) = 0.586（`|Δ| = 0.001`）、σ(B) = 0.375，签名成立，细化量子后**更干净**。
- **结论修正（U2 类型 E）**：推论 4.4 → **4.4′**。最优注入的主要收益**不是缩短收敛时间，
  而是降低伪收敛率**（每个预算档位 F2 都系统性更低：5.3 vs 8.2、6.0 vs 7.8、4.8 vs 6.8）；
  收敛速度优势只在资源紧张时出现（B=10 处 +9.7 pp）。
- **预言撤回**：R18-b（内约束强度的剂量-反应）**撤回**。拆解后发现 κ 的两个子机制都不成立——
  停止标准那半是**同义反复**（τ 低于接受阈值当然不收敛），注入策略那半在主参数下**完全无效应**
  （68.0 / 70.0 / 68.3 / 70.3 / 69.8，无差异）。**内约束强度尚无有效操作化**，登记为缺口 **B-B21**。

#### fix: 主实验 seed 派生不可复现（入库前实测发现）

- 原 `run_arm` 用 `seed + hash(arm) % 9973` 派生各臂种子。Python 的 `str` 哈希**受 `PYTHONHASHSEED`
  随机化影响**，同一 seed 在不同进程产出不同结果——入库前实测：连续两次运行 `r18_dual_constraint.py`
  得到的 `r18_results.json` **md5 不一致**。
- 改为 `zlib.crc32(arm.encode("utf-8"))`，修复后两次运行 md5 一致。
- **本目录下所有 json 均为修复后重跑的产物**，与初版报告数字有小幅差异（如 base 真收敛 69.6 → 71.1），
  **结论不变**。已用 `test_seed_derivation_is_process_stable` 与
  `test_source_does_not_use_builtin_str_hash` 锁死，防止回退。

#### fix: 更正初版报告中的一处引用错误（勿删此条记录）

初版报告写"高分辨率复现 base σ = A σ = 0.581 精确相等"。该引用**是错的**——
0.581 取自 `r18_check.json` 的 κ 扫描（`N_FEAT=16, κ=1.0` 那一行的 σ），**不是 base/A 的对照值**；
`r18_sweep.json` 的鲁棒性网格虽含 `N_MACRO=40`，但 `N_FEAT` 仍为 8，**当时根本没做过高分辨率的 base/A 对照**。
新增 `r18_highres.py` 补做该对照（真实值 0.587 / 0.586）。结论方向未变，**但引用链此前是断的**。

同理，"σ(A) 与 σ(base) 精确相等"的措辞收敛为**近似相等**（`|Δσ| ≤ 0.006`）——
当时的"精确相等"是特定 seed 下的巧合，**不是结构性质**。

#### test: `tests/test_r18_dual_constraint.py`（7 项）

只锁**判别签名与可复现性**，不锁具体数值（数值随参数变化，锁死会阻碍任何参数探索）：

- `test_outer_failure_is_F2_dominant` / `test_inner_failure_is_F1_dominant`——两类失效分工
- `test_signature_sigma_outer_unchanged_inner_dropped`——**签名 R18-a**，且要求内外分离度相差一个数量级
- `test_F1_F2_are_derived_not_assigned`——F1/F2 由 σ 路径导出；三态在操作化层面均可达
- `test_seed_derivation_is_process_stable` / `test_source_does_not_use_builtin_str_hash`——可复现性防回退
- `test_signature_holds_at_higher_resolution`——排除"签名只是 σ 离散化伪像"

#### 已知局限（新增，详见 README「已知局限」6–8 条）

- R18 是 **Python 层面的操作化**，与 `formal/` 的 Lean 4 定理之间尚无对应（**B-B17**）。
- 宏库是人造的（原型 + 扰动），不是自然语言；结论关于**机制**，不关于语言。
- 单输入模型：外失效时系统无交叉检验来源，F2 率是**上界**。
- B 臂 σ 标准差 ≈ 0（恒定 0.375），动力学被停止规则压扁，"F1 = 100%"比看起来更依赖建模选择。

#### 理论侧缺口（详见 `docs/xs-framework/README.md` 第四节）

未闭合（高）：**B-B19**（信息层 U 与动力学层 Ψ 的桥未证，在此之前两层不得互引）、**B-B11**。
未闭合（中）：B-B16 / B-B17 / B-B21。已闭合：B-B6 / B-B10（主要部分）/ B-B13 / B-B14。

## v1.6.1 — 2026-10-18
### R3-v2 案例库定量升级（S3 流感 / S5 鳕鱼，修订登记 §17.7）
- feat: `src/r3_case_upgrade.py` — 公开来源数据集（CDC 逐季 VE 2004-2026、WHO H3N2 组分
  2004-2026、WHO 2014-15 错配时间线、DFO 北方鳕 SSB/死亡率锚点）+ 七条预注册判据
  （C-S3-1..3 / C-S5-1..4）+ 判定函数
- feat: `experiments/sim_r3_case_upgrade.py`；图 `figs/r3_case_upgrade.png`
- test: `tests/test_r3_case_upgrade.py`（14 项，含数据完整性防回写漂移；全套 134 项通过）
- findings（S3）：漂移季 2014-15 VE=19% 为紧 CI 季系列最低（z=−2.76，均值 0.45 倍）；
  WHO 换株率 15/21=71%（瞬时读法为现役协议）；2014-15 窗口逐点成立
  （决策时点 >90% 匹配 → 季中 <50% 抑制 → VE 19%，当季 9 月改荐 Switzerland 株纠正）。
- findings（S5）：耗竭期 F>M（1950s-80s）+ 崩溃 5 年 −95%；禁捕后低于 LRP 24 年
  （吸收态实测）；迟滞比 5.8；崩溃后 M 均值 0.47=2.35× 崩溃前假设、Capelin 仅 10%
  （再生机制受损 = S5 的 phi_min 通道的现实对应）。

## v1.6.0 — 2026-10-18
### 五项遗留动态事项一次清账（T1-d / 方案甲 / R3-v2 / R4-LE / ΔT 操作化）
全部经「预注册判据 → 模块 → 实验脚本 → 单元测试 → 图」流水线完成；全套回归 120 项通过。

#### T1-d 表征分化仿真（第八节 R8 遗留）
- feat: `src/t1d_representational_diff.py` — 双阱势 + 多单元群模型（300 单元、涨落驱动
  η_eff=σ√(1+J)、旁观维度弱噪声）；五判据边界诊断（BC/sep/abs_sep/straddles/symmetry）
- feat: `experiments/sim_t1d_representational.py`；图 `figs/r8_t1d_representational.png`
- test: `tests/test_t1d_representational.py`（12 项）
- findings: 三子命题全成立——翻越率=0.0110+0.0086·J（R²=0.970，时点归流）；
  a=0 阴性对照任意 J 不生成边界（存在归势）；f₊∈[0.49,0.51]（方向归流）。
  **新发现**：BC 随 J 从 0.739 降至 0.522——通量过大边界不可维持（A11b 非单调性证据）。
- 建模教训（已入模块文档串）：v1 乘性可分化性谱只改方差不改均值结构；
  v2 单轨迹不锁阱；v3 多单元群 + 涨落驱动才成立。

#### 方案甲（A0 回归收敛）仿真（R7 遗留）
- feat: `src/a0_bootstrap_convergence.py` — 层级回归链 ε_k=g·ε_{k−1}+e_k、
  C-1/C-2/C-3 三档判据、不动点迭代器、R7 四类人工系统算子模型
- feat: `experiments/sim_a0_bootstrap.py`；图 `figs/r7_plan_a_convergence.png`
- test: `tests/test_a0_bootstrap.py`（9 项；含半正态口径修正——tail_mean 对照
  E|ε|=σ·√(2/π)/√(1−g²)，二阶矩对照方差）
- findings: 不动点迭代观测收敛比精确等于 g；四类系统定位——类1 保距（C-2）、
  类2 弱收缩到常驻偏置、类3 扩张（C-3，方案甲被证伪）、类4 强收缩但不动点是
  适应峰≠A0（收敛≠回归 A0，支持 R7「丙为主乙为辅」）。

#### R3-v2：C1+C2 充分性补案（R3 遗留）
- feat: `src/r3_sufficiency_probe.py` — 预注册判据 + 五类构造
  （S1 累积损伤 / S2 竞争置换 / S3 时序错配 / S4 中介污染 / S5 势耗竭）
- feat: `experiments/sim_r3_sufficiency.py`；图 `figs/r3_sufficiency_probe.png`
- test: `tests/test_r3_sufficiency.py`（11 项，含参数敏感性与不合法构造识别）
- findings: **充分性被证伪**——五类「单次全合规、序列全失败」构造全部判强反例。
  C1+C2 作用域须限定为**单次交互的准入判据**；序列成败须追加序列条件。
  **条款修正**：S3 逼出 C2 必须按**瞬时读法**（|vL−vE| 对照绝对基准）判定——
  增量读法下漂移可在全程合规中越界。

#### R4 分量归类的 LE 内核验证（R4 遗留）
- feat: `src/r4_component_attribution.py` — F2 守恒残差采集（多初值集成破共线性 +
  C1 失败步屏蔽）、分解回归（α=−(ratio−1), β=−base）、K1 状态相关性 / K2 修复检验 /
  K3 滤波检验（δ 翻号符号不变性）、三域归类
- feat: `experiments/sim_r4_attribution.py`；图 `figs/r4_component_attribution.png`
- test: `tests/test_r4_attribution.py`（12 项，含纯失配/纯底噪对照与反向判别力）
- findings: 复合残差无歧义分解，α、β **精确**恢复解析真值（R²=1.000）；
  失配分量 |ρ|≈0.95–0.98 → 接口域（修复后消除率 100%）；底噪分量符号不变 →
  载体域（滤波消除）。**可辨识性前提**：分解要求 δ 充分变号且屏蔽 F1 步——
  这是判据的适用条件，须写入条款。

#### ΔT（时间势差）/ TP 谓词操作化（V3.3.1 §G4 登记的缺口）
- feat: `src/tp_temporal_potential.py` — ΔT=|ω_L−ω_E|/(ω_L+ω_E)；耦合相位振子
  （Adler 方程）；观测量 = 相对相位漂移率；P1–P5 预注册判据
- feat: `experiments/sim_tp_temporal_potential.py`；图 `figs/g4_tp_temporal_potential.png`
- test: `tests/test_tp_temporal_potential.py`（10 项）
- findings: ΔT→0 通道通过——解锁区漂移率与 Adler 精确解 √(Δω²−4K²) 一致
  （误差 < 3%）；|Δω|≤2K 锁定（势差被接口吸收）；K=0 阴性对照漂移恒=Δω
  （锁定是耦合诱导的）；ΔT 与 d 独立可变（P4，TP 非冗余谓词）。
  框架语义：锁定 = 时间势差在接口内被消化；漂移 = 势差外显为不可逆时序流
  （A11a/b 势-流分离的时间版）。

## v1.5.0 — 2026-10-18
### G1-a：LE 公理体系的 Lean4 形式化工程（formal/）
- feat: `formal/Formal/Basic.lean` — 类型/谓词骨架（Sys/Scale/Tau、leTau 预序、约 40 个谓词签名）
- feat: `formal/Formal/Axioms.lean` — LE 公理 L1 层编码：SameScale 等价关系性、耦合充要条件
  `coupling_iff : Coupled l e ↔ SameScale l e`、不匹配仅单向扰动、M1/M2、F1–F16、E1–E7、H_LE/H_EL
- feat: `formal/Formal/Corollaries.lean` — 可证推论 9 条（tactic 证明，禁 sorry 纪律），含
  mismatch_forbids_coupling（iff 逆否）、aging_mechanism_exists、life_produces_waste_heat、
  coevolution_bidir、life_writes_back 等
- feat: `formal/Formal/AttackSurface.lean` — N1 否定式攻击面分析 + 一致性守卫
- docs: `formal/GAPS.md` — 编码暴露的公理缺口 B1（衰老缺"阈值必破"桥接）/B2（演化桥接，
  须带资格条件字段 C_q，与 R5 逐语境结论的兼容性待裁决）/B3（吸引子分类超出一阶骨架）/
  C4-full（耗散结构完整版缺合系统孤立性）；候选桥接公理+相对一致性注记
- findings（不依赖编译结果）：归档版第 7 节"推论全量可演绎导出"声明**被证伪**——
  6 条声明推论中 3 条依赖归档版未登记的隐含前提
- ✅ **类型检查已通过**（2026-10-18）：四文件全部经 Lean 4.14.0 内核编译通过（linux_wasm32 官方
  构建 + Node 运行时，native lake 被本机程序黑名单拦截；同一 elaborator，结论等价）。零错误、
  零 sorry/admit。核验中修复两处真实错误：模块文档（/-! -/）挡 import（复合命令序违例）×2、
  evolution_materials 漏 intro 第三假设 ×1。本工程现可作为 N1 否定式的公开攻击面引用。

## v1.1.0 — 2026-10-17

中继链（命题 2′ 中继部分）落地；两项理论修正由仿真裁决。

### Added
- **feat**: `src/relay_chain.py` — 串联中继链模块：`RelayChain`（多段输运/自演化）、
  `simulate_transport`（夹持端点稳态传输率测量）、`min_relay_count`（桥接给定端距的最小中继数）、
  `leg_overlap`（与内核同式的重合度）
- **feat**: `experiments/sim_relay_recovery.py` — 三联图：F1 恢复对比 / 复合律裁决 / 失效模式；
  输出 `figs/relay_recovery_and_law.png`
- **test**: `tests/test_relay_chain.py` — 10 项测试（复合律、F1 阻断、F2 污染的解析值核对、
  中继数下限、恢复成功/失败）
- **docs**: `docs/experiment_manual.md` §6「中继链实验」

### Changed（理论修正——由本版仿真裁决）
- **命题 2′ 复合律修正**：串联中继的传输强度服从**调和/短板律** T = 1/Σ(1/dᵢ)（实测逐点吻合，
  n=2 时 0.45 vs 乘积律 0.81）；∏(pᵢqᵢ) 仅适用于"各段独立必要的并联/事件语义"（路径存在概率）。
- **C2 角色修正**：F2 段不阻断主模态传输（T 不变），只注入解析可预测的 ghost 污染
  （0.25/(1+1.6)，随链长增大）→ C2 是**可信度**条件，不是**连通性**条件。
- **新预测**：中继数下限由弱端尺度决定（gap 1.09 需 3 段；gap 0.6 需 1 段）——I1 可构造性的
  定量代价，与 MMS 分级桥接（MAAD/bridging domain）同构。

## v1.0.0 — 2026-10-17

LE（Life-Environment，生命–环境同尺度耦合）原型首次落代码。此前仓库仅有 README.md。

### Added
- **feat**: `src/le_coupling_kernel.py` — LE 双投影子系统耦合动力学内核
  - `calc_d_overlap`：边界重合度标量 d（无量纲，高斯核，双侧尺度归一化）
  - `judge_interaction_conditions`：C1/C2 双条件判定 → high / middle / low 三区 + F1/F2 失效标记
  - `compute_coupling_flux`：载体通量；F1 置零、正常耦合作用-反作用守恒、F2 副模态幅值失配产生 ghost 伪力
  - `LEProjectionSystem.test_D_encapsulation`：A10 边界封装检验（滚动窗口 + 双峰系数 BC + 偏度 + 方差界）
  - 势-流分离：`semantic_potential_grad` 只进耦合增益，不进入状态导数（A11-R）
- **feat**: `experiments/sim_le_coupling.py` — 四联时序图（L/E 状态、d、Var(D)+ghost、三区/F1F2 标记），自动保存 `figs/le_coupling_timeseries.png`，附分段统计摘要
- **test**: `tests/test_le_kernel.py` — 19 项单元测试：d 值域与单调性、C1/C2 三区、D 封装（单峰通过 / 多峰与大方差拒绝 / 样本不足不误判）、F1 零通量、F2 伪力与守恒破坏、C2 兼容时通量守恒、历史登记完整性、势不进入状态方程、耦合确实改变轨迹、种子可复现、属性别名
- **docs**: `docs/experiment_manual.md` §5「LE 生命–环境耦合仿真实验」：参数表、分段逻辑、结果解读、理论对照
- **docs**: README 补全此前缺失的 `## Repository structure` 章节与仿真脚本清单

### Fixed（相对初版方案的设计错误）
- 修正 **阶段顺序**：原方案「先 F1 30–45s，再 F2 45–60s」在动力学上不可能产出 F2 数据——一旦跌入 F1，`d` 已低于 C1 阈值，之后关闭 C2 只会继续被判为 F1。改为 S1 稳定 → S2 F2 → S3 F1 → S4 恢复尝试
- 修正 **状态采集**：不再在绘图段复制内核状态，`history` 直接登记 `L_state`/`E_state` 快照
- 修正 **属性名不一致**：`semantic_potential_grad` 与 `semantic_pot_grad` 双名并存（后者为 property 别名），避免 `AttributeError`
- 修正 **D 指标口径**：封装检验改为滚动窗口（默认 200 步），避免「历史混合分布」把早期稳定态与后期失效态混成一个假多峰
- 修正 **启动伪影**：样本数不足时返回 `status="warmup"`，不计入 M-不完备

### Known limitations
- 参数未标定到真实生态数据；3 维状态是占位实现，不是 LE 体系 F1–F17/E1–E7 的完整映射
- F2 的 ghost 伪力是**现象学建模**（幅值失配 + 底噪项），尚未从界面算子（提升/限制算子不相容）解析推导
- F1 的「吸收态」性质部分来自 `d` 乘在通量里这一建模选择，需与真实生殖隔离的不可逆性做对照检验，不能直接当作理论预测
- 认知架构（流形投影 + 非对称势阱，T1–T6）的仿真代码仍不在本仓库

## v1.2.0 — 2026-10-17
### R2 跨域 D 形态审计
- feat: src/d_form_audit.py — D 形态统计（绝对/申报双标度）、三参数合成器、秩相关
- feat: experiments/r2_constraint_ordering.py — 物理/生物/人文三域对照；表 1 排序、表 2 等峰度线、表 3 可辨识性；图 r2_d_form_ordering.png
- test: tests/test_d_form_audit.py（8 项；全套 37 项通过）
- findings: 核心宽度排序与 τ 一致但属构造性；峰度排序与 τ 不一致且非单调；等峰度线对角 → 跨域不可辨识；σ 倍数不是跨域可比单位；C2 获操作化定义（跨实现一致率）

## v1.3.0 — 2026-10-17
### R6 A10 可观测现象清单（预注册 + 判别力检验）
- feat: src/a10_observables.py — 预注册判据表（P2–P5，绝对标度）、粒度扫描（残差–粒度标度律）、持续性检验（3 窗口 ≥2）、三联排除判定（区域门 → 修接口 → 标度律）
- feat: experiments/r6_a10_observables.py — 六场景判别力测试（OK/F2/F1/SUB_OU/SUB_1F/BIAS）+ τ 扫描判别力边界；图 r6_a10_observables.png
- test: tests/test_a10_observables.py — 16 项（含判据表冻结快照、白噪声 β≈−0.5、1/f β≈0、顺序置换、位置不变盲区；全套 53 项通过）
- findings: 五个可判场景全部裁决正确，A10-03 解除；D 形态判据无判别力（三失效同形），三联排除为必要结构；排除顺序 load-bearing（F2 反序误判）；P5 是"特征尺度是否在扫描窗下部"探测器（τ≈25–128 步的有限尺度被误判为坍塌）；盲区 G2 证实——系统性偏置须靠跨实现一致率（R2 的 C2）捕获

## v1.4.0 — 2026-10-17
### R5 推论 3 的 A9 资格审计
- feat: src/a9_qualification.py — Wright–Fisher 仿真核（扩散解析式锚定）、d 交叠构念、判据②确定性/统计/干预三读法、判据③冻结迁移、判据①闭合钳制、资格裁决
- feat: experiments/r5_a9_qualification.py — Ns 梯度四语境审计 + 迁移泄漏 + 序层方向（含 Ns→0 极限点）；图 r5_a9_qualification.png
- test: tests/test_a9_qualification.py — 13 项（P_fix 解析式锚定、中性基线、干预失敏带、闭合泄漏、裁决结构；全套 66 项通过）
- findings: 原降级理由（确定性不可预测）不成立——统计版判据②全境 d≥0.93；真正边界是干预失敏（Ns≲5 时 eff<0.2）与冻结参数跨数量级迁移崩塌（③）；正确表述为逐语境宏 M「Ns≳5 ∧ m=0 ∧ (s,N) 申报」；"无向"是 Ns→0 极限性质（D_seq: 0.01→0.90）；A9 需增补资格条件字段 C_q
