# X-S-MERA

Manifold Projection and Asymmetric Potential Well: Cognitive Dynamics Framework

> This repository contains the reference simulation code for the X‑S‑MERA theoretical framework.
> Author: [zjg]
>
> Current version: **v1.6.40**（2026-10-29） · Tests: **393 passed** · License: **MIT**
>
> 测试数为**含 numpy 环境**（CI 的 `requirements.txt`）下的全套计数；本机若缺
> numpy，11 个历史模块（LE 内核 / relay / R2–R4 / A0 / A9 / A10 / T1-d / ΔT）
> 会在 collection 阶段报 `ModuleNotFoundError`（**环境问题，非回归**）。
> 另有 R29 单跑约 78 s，合并跑需注意超时。

## Abstract

X‑S‑MERA interprets AI hallucination, reasoning oscillation and cognitive stagnation as **endogenous
information residuals caused by manifold dimensionality‑reduction projection**, rather than purely
external engineering defects. The LE (Life–Environment) prototype extends the same formal machinery
(边界重合度 d、D 分布封装、C1/C2 双条件、F1/F2 失效) to coupled subsystem dynamics.

## Repository state（重要，诚实声明）

| 内容 | 状态 |
| --- | --- |
| LE 生命–环境耦合内核 `src/le_coupling_kernel.py` | ✅ 已实现，19 项单元测试全过 |
| LE 耦合仿真 `experiments/sim_le_coupling.py` | ✅ 已实现，出图见 `figs/` |
| 中继链 / R2 / R6 / R5 审计模块 | ✅ 已实现（v1.1.0–v1.4.0） |
| 五项遗留动态事项（T1-d / 方案甲 / R3-v2 / R4-LE / ΔT 操作化） | ✅ 已实现（v1.6.0） |
| R3-v2 案例库定量升级（S3 流感 / S5 鳕鱼，公开数据集 + 七条预注册判据） | ✅ 已实现（v1.6.1），14 项测试；**全套回归 134 项通过** |
| 认知架构仿真（流形投影 + 非对称势阱，T1–T6 自证伪测试） | ✅ **已落地（v1.6.40）**：`src/kernel.py`（见下行）；内核级自证伪 T-affine/T-gauge 已实现，T1–T6 实验臂登记下一迭代 |
| **认知架构内核 `src/kernel.py`** | ✅ 已实现（v1.6.40）。状态(E,U)+规范固定参数+仿射积分器+A_D⁺的g(U)结构约束；T-affine/T-gauge 内核级自证伪；D-SCALE 三值判别接口（R34 禁区 {|β̂|≤z·SE}）；R23 输出适配器（to_r23_json + ρ_res>ρ_exp 走向）。14 项测试（`tests/test_kernel.py`）。诚实边界：数值原型非经验验证；A_D 未证（B-B22）+单一驱动附款（S2，调用方契约）；快思维本体不可测（Q0）；g(U)形式是约定（R45）；读数条件于 A_D 与声明的尺度 |
| LE 公理体系 Lean4 形式化工程 `formal/` | ✅ 已通过 Lean 4.14.0 内核类型检查（wasm 运行时；native lake 见 formal/README.md） |
| LICENSE | ✅ MIT（见 `LICENSE`）；release 待发布 |
| CITATION | ✅ `CITATION.cff` 已补 |
| **XS 框架理论文档集** `docs/xs-framework/`（v3.1 认知框架主线 + v3.0 统一版 + 8 月稿与审计） | ✅ 已入库（v1.6.6），**19 篇**；索引见该目录 `README.md` |
| **R18 双约束失效可分离性实验** `experiments/r18_dual_constraint/` | ✅ 已实现（v1.6.2），7 项测试；报告见该目录 `report.md` |
| **R19 U/Psi 桥一致性检验** `experiments/r19_bridge/` | ✅ 已实现（v1.6.3），7 项测试 |
| **R21 内约束强度 = O_adm 剪枝量** `experiments/r21_constraint_specificity/` | ✅ 已实现（v1.6.3），7 项测试 |
| **R27 双相形状：X 单调 / M 单峰** `experiments/r27_biphasic/` | ✅ 已实现（v1.6.4），14 项测试。峰位命中解析预测、**峰序吻合率 1.000**；**单峰性须粗粒化后检验** |
| **R28 微观回落 / 多峰 / 局部-全局收敛** `experiments/r28_micro_reversal/` | ✅ 已实现（v1.6.5），15 项测试。回落率解析律命中；**路径层单峰率 0.000（本体性，非测量误差）**；局部收敛 100% vs 全局收敛 13–34% |
| **R29 峰距标定程序** `experiments/r29_peak_distance_calibration/` | ✅ 已实现（v1.6.6），15 项测试。公式命中（误差 ≤1 步=离散化界）；估计量无偏；**异 K 必须归一化**（raw 差因子 2.000）；**共享噪声下峰序零判反**；U−E 恒等式残差 5.7e-14 |
| **R-042 D-SCALE 参数敏感性扫描** `experiments/r42_dscale_eta/` | ✅ 已实现（v1.6.14）。三类场景 × η∈[0.005,1.0] 判定读数**正确区分且对 η 稳定**，g 不进入判定流程 ⟹ **判定程序层已过**（g、τ\*/η 均未越界，U5 通过）。**扫描副产品 R-043 曾疑"τ\* 落带区间几乎为空" ⟹ 已解决（v1.6.23）：加密 η 网格（200 点）后落带区间非空，原疑为网格假象**；η 取值原则由 **R-044（v1.6.32）** 锚定 ⟹ **B-B44 闭合**（见下两条） |
| **A5 generics 一般性证明** `experiments/a5_generics/` | ✅ 已实现（v1.6.15）。四条分类规则（楔形/峰序/四级可辨识/共线反解）在 A_D 原型内满足函数性 genericity，720 实例零违反；代数证明 + 原型枚举，**非 Lean、非现实断言**（见 `proof.md`） |
| **CI 强制门** `.github/workflows/xs-gates.yml` | ✅ 已交付（v1.6.14）且**已在 GitHub 启用**（2026-10-27 确认）。`xs_trace_gate.py`(C1.7/C1.8) + `xs_formal_checks.py`(C1/C2) + `xs_boundary_map_check.py`（A10 残留，v1.6.38） + `pytest tests/` 四道合并门槛 |
| **A7 核心最小化** `experiments/a7_core_min/` | ✅ 已实现（v1.6.16）。47 记号→11 核心（根 10+约定输入 1），36 导出全部可达，压缩率 76.6%；依赖结构层面，非概念独立性 |
| **R18 的 Lean 命题规格** `formal/Formal/DualConstraint.lean` | ✅ 已接入（v1.6.18）。认知双约束/R18 Horn skeleton 编码，lake build 编译通过、无 sorry；见 `formal/GAPS.md` 条目 C1 |
| **A13 第一例对应探索** experiments/a13_first_case/ | ✅ 已实现（v1.6.17）。换工作决策：D-DRILL 钻研/空洞相位 + D-CONV 收敛/发散方向 + D-IDENT 分级，对照无框架程度描述；演示性映射，非现实测量 |
| **B 组方法准备**（B3/B1/B2，外部数据前可执行部分） | ✅ B3 残差反演 g（v1.6.19，4 种 g corr=1.0）；B1 四臂闭合（v1.6.20，自洽 1e-14）；B2 R31 守卫（v1.6.21，正常 0.000 / A_D⁺ 违反 1.000）。**非经验验证**；R43′（独立 Ψ）、真人测量、跨人群仍属 B 组待外部数据 |
| **v1.6.14–1.6.21 新增实验的 pytest 覆盖** `tests/test_{a5,a7,a13,b1,b2,b3,r42}_*.py` | ✅ 已补齐（本地未提交）。7 套共 **67 项**结构性测试（锁结构预言、不锁数字，仅标准库，1.2 s 跑完）；**CI 强制门 `xs-gates.yml` 的 pytest 步骤自动覆盖**，无需改 workflow |
| **R8 / R11 / R44 的 pytest 覆盖** `tests/test_{r8_period_doubling,r11_pseudo_drift,r44_eta_principle}.py` | ✅ 已补齐（v1.6.34）。3 套共 **15 项**（R8 8 项锁 M1 周期倍增 V1–V4；R11 5 项锁推论 4.1 的 legal/violate 双臂区分；R44 2 项锁 τ\* 对 η 单调 + 可靠子区间非空）。均**调底层函数复现结构预言、不调 `main()`**，零文件写入副作用；仅标准库 |
| **R8-B9 唯一性归约** `experiments/r8_gain_stepchain/r8b_uniqueness.py` | ✅ 已实现（v1.6.37）。机制消融四判据 U1–U4 全真 ⟹ **框架内唯一性成立**：A11-R 下唯一能改变混沌拓扑的势机制是 M1；M2 合法形式归约 M1（增益不对称=多候选有效步长调制）、非法形式（状态偏向）被排除；M3 强制归宿只压缩混沌寿命（300000→51）、不移动 g_chaos（恒 3.57）⟹ **B-B9 推到位移**（框架内唯一性已证；跨框架全局唯一性转外部开放）。诚实边界：S2 最小离散原型，非全局唯一性主张 |
| **边界地图机检** `src/xs_boundary_map_check.py` | ✅ 已实现（v1.6.38）。A10 残留落地：`boundary_map.json` 从"仅登记"推进到"被机器看住"——逐版本结构/枚举/解决指向校验，已接入 CI 门第四步；`boundary_map.json` 升级 v0.6；11 项结构测试（`tests/test_xs_boundary_map_check.py`） |
| **B-B32 震荡判据一致性** `experiments/b32_shape_prob_consistency/` | ✅ 已实现（v1.6.38）。S2 合成检验推论 6.11 形状判据 vs 主文 U 概率判据，**一致率 0.75、双向均有反例，一致性不成立**：反复重开震荡中 U 被重开抬高不满足"U 仍在下降"（概率判据漏判）；从未到 m=1 的未决振荡中无局部收敛可重开（形状判据漏判）⟹ **二者互补而非冗余，需合并覆盖**（08 supp7 缺口表 B-B32 → 已检验）。诚实边界：判据操作化为建模选择 |
| **礁石 R-043（τ\* 落点）** `docs/xs-framework/reef_shturl.json` | ✅ **已解决（v1.6.23，`experiments/r43_tau_star/`）**。原疑：τ\* 落在 (τ_fine,τ_coarse) 内的 η 区间在该网格上几乎为空。**加密 η 网格（200 点）后三类场景在窗口 (8,128) 内均有非空落带区间**（A[0.031,0.046]／B[0.031,0.071]／C[0.061,0.086]）——此前"几乎为空"是**网格假象**。η 取值**场景相关**（三区间交集为空），不由单一全局值给出 |
| **R-044 η 取值原则锚定** `experiments/r44_eta_principle/` | ✅ 已实现（v1.6.32）。落带区间内 Bootstrap 压出**可靠子区间**（落带率≥90%：A[0.031,0.046]／B[0.032,0.071]／C[0.061,0.086]）；**取值原则 P1–P5 定稿**（`33-v3.4-eta-scaling-principle.md`）⟹ **B-B44 闭合**。诚实边界：可靠子区间为**原型内**结果，REL_TH=0.90 为约定水平；η 场景相关，**不存在全局单值**（须按场景显式声明） |
## Repository structure

```
x-s-mera/
├── src/                          # 唯一真值源（experiments/ 与 tests/ 只依赖此处公开接口）
│   ├── le_coupling_kernel.py     # LE 生命-环境耦合动力学内核（d/D、C1+C2、F1/F2、势-流分离）
│   ├── relay_chain.py            # 中继链（命题 2′ 中继部分；调和/短板律）
│   ├── d_form_audit.py           # R2：跨域 D 形态审计
│   ├── a10_observables.py        # R6：A10 可观测现象清单（三联排除）
│   ├── a9_qualification.py       # R5：A9 资格条件审计（Wright–Fisher 核）
│   ├── t1d_representational_diff.py  # T1-d：表征分化（双阱势+多单元群）
│   ├── a0_bootstrap_convergence.py   # 方案甲：A0 回归收敛（C-1/2/3 判据）
│   ├── r3_sufficiency_probe.py       # R3-v2：C1+C2 充分性反例全谱
│   ├── r3_case_upgrade.py            # R3-v2 案例升级：S3/S5 公开数据集 + 预注册判据
│   ├── r4_component_attribution.py   # R4：A8/D17 分量归类（K1→K2→K3）
│   ├── tp_temporal_potential.py      # ΔT：时间势差 TP 谓词操作化（Adler 锚定）
│   └── kernel.py                 # ⏳ 计划中：认知架构内核（快慢耦合非对称势阱），当前不在仓库
├── experiments/                  # 每个 src 模块对应的出图脚本（sim_*/r*_*.py）
│   ├── a5_generics/              # A5：generics 一般性证明（1 脚本 + 结果 json + proof.md）
│   ├── a7_core_min/          # A7：核心最小化（1 脚本 + core_table.json + README）│   ├── r18_dual_constraint/      # R18：双约束失效可分离性（5 脚本 + 5 json + report.md）
│   ├── a13_first_case/       # A13：第一例对应探索（1 脚本 + 对照 json + README）│   ├── r19_bridge/               # R19：U/Psi 桥的双速率分解一致性检验（1 脚本 + json + report.md）
│   ├── r21_constraint_specificity/  # R21：内约束强度 = O_adm 剪枝量（3 脚本 + 3 json + report.md）
│   ├── r27_biphasic/             # R27：双相形状 X 单调 / M 单峰（1 脚本 + json + report.md）
│   ├── r28_micro_reversal/       # R28：微观回落 / 多峰 / 局部-全局收敛（1 脚本 + json + report.md）
│   ├── r29_peak_distance_calibration/  # R29：峰距标定程序（1 脚本 + json + report.md）
│   ├── r30_gauge_fixing/         # R30：A_D 规范固定恒等式（1 脚本 + 结果 json）
│   ├── r33_clevel_closure/       # R33：C 级闭合实测（1 脚本 + 结果 json）
│   ├── r34_vertex_band/          # R34：ρ=1 顶点容差带（1 脚本 + 结果 json）
│   └── r42_dscale_eta/           # R-042：D-SCALE 对 η 参数敏感性扫描（1 脚本 + 结果 json + 运行留档）
├── tests/                        # 301 项单元测试，全过（23 个测试文件，v1.6.14 实测）
├── docs/
│   ├── experiment_manual.md      # 实验手册（§5 LE 耦合实验、§6 中继链实验）
│   ├── xs-framework/             # XS 框架理论文档集（16 篇，含 v3.1 认知框架主线）
│   │   └── README.md             # 索引：版本线、U1 取代纪律、缺口表、与代码层对应
│   └── summary_theory.md         # ⏳ 计划中：理论摘要，当前不在仓库
├── formal/                       # G1-a：LE 公理体系 Lean4 形式化（L1 对象层 Horn 骨架编码）
│   ├── Formal/                   # Basic（谓词骨架）/ Axioms / Corollaries / AttackSurface
│   ├── GAPS.md                   # 编码暴露的公理缺口 B1/B2/B3/C4-full 与候选桥接
│   ├── lakefile.toml / lean-toolchain
│   └── README.md                 # 构建命令与状态声明（✅ 类型检查已通过）
├── figs/
├── requirements.txt              # numpy>=2.0 / matplotlib>=3.8 / pytest>=8.0
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
└── README.md
```

## Quickstart

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt        # Windows
# source .venv/bin/activate && pip install -r requirements.txt   # Linux/macOS

python experiments/sim_le_coupling.py     # 出图 + 分段统计摘要
python -m pytest tests/ -v                # 393 tests（R29 单跑 ~35 s，全套 ~70 s）
python experiments/r18_dual_constraint/r18_dual_constraint.py   # R18 主实验（仅标准库）
python experiments/r19_bridge/r19_bridge.py                     # R19（仅标准库）
python experiments/r21_constraint_specificity/r21_constraint_specificity.py  # R21（仅标准库）
python experiments/r27_biphasic/r27_biphasic.py                              # R27（仅标准库）
python experiments/r28_micro_reversal/r28_micro_reversal.py                    # R28（仅标准库）
python experiments/r29_peak_distance_calibration/r29_peak_distance_calibration.py  # R29（仅标准库）
```

复现须知：

- 需要 Python ≥ 3.10（依赖 `numpy>=2.0`）。
- `tests/test_r3_case_upgrade.py` 含数据完整性校验（防回写漂移），耗时明显高于其余文件；
  只做快速自检可先跑 `python -m pytest tests/ -v -x --ignore=tests/test_r3_case_upgrade.py`。
- 全套 **301 项**（v1.6.14 实测 `pytest tests/` → 301 passed；含 v1.6.6 之后新增的 r30_gauge_fixing、r33_clevel_closure、r34_vertex_band 及校验器测试）。原 141 项分发（v1.6.6）：
  le_kernel 19 / a10 16 / a9 13 / r3_case_upgrade 14 / r4 12 /
  t1d 12 / r3_sufficiency 11 / relay_chain 10 / tp 10 / a0_bootstrap 9 / d_form_audit 8 /
  r18_dual_constraint 7 / r19_bridge 7 / r21_constraint_specificity 7 /
  r27_biphasic 14 / r28_micro_reversal 15 / r29_peak_distance_calibration 15。
- `experiments/r18_dual_constraint/` 下的脚本**只依赖标准库**，不需要 numpy；
  运行时把 json 写进当前工作目录，故请在该目录内执行（见 `report.md` §7）。
- **CI 强制门（v1.6.14 · A10）**：`.github/workflows/xs-gates.yml` 把 `xs_trace_gate.py`（C1.7/C1.8）、
  `xs_formal_checks.py`（C1/C2）、`pytest tests/` 串成三道合并门槛，任一 error 即阻断。
  本地等价验证：`python src/xs_trace_gate.py && python src/xs_formal_checks.py && python -m pytest tests/ -q`，须全退出码 0。

## 理论 → 代码对照

| 理论条款 | 代码实现 | 检验 |
| --- | --- | --- |
| d：边界重合度标量（无量纲，非时空尺度比） | `calc_d_overlap` | `test_d_overlap_calc_range_and_monotonicity` |
| 命题 2′：C1 重合度 + C2 接口兼容，三区 | `judge_interaction_conditions` | `test_C1_C2_judge_three_region` |
| F1：无公共细化 → 直接交互禁止，通量置零 | `compute_coupling_flux` | `test_F1_branch_flux_zero` |
| F2：有重合、接口不兼容 → ghost 伪力（礁石伪收敛） | `compute_coupling_flux` F2 分支 | `test_F2_branch_ghost_flux` |
| A10：D 可封装 ⇔ 边界宏；偏离 ⇔ M-不完备 | `LEProjectionSystem.test_D_encapsulation` | `test_D_encapsulation_*` |
| A11-R：势只调制增益，不进入状态方程 | 通量公式中势仅出现在 gain | `test_potential_never_enters_state_equation` |
| **双约束两类失效分工**：外约束失效（对象被替换）→ F2；内约束失效（任务被替换）→ F1 | `experiments/r18_dual_constraint/` 四臂对照 | `test_outer_failure_is_F2_dominant` / `test_inner_failure_is_F1_dominant` |
| **判别签名 R18-a**：σ(外失效) ≈ σ(base) 且 σ(内失效) < σ(base) | 同上，σ 为中介变量 | `test_signature_sigma_outer_unchanged_inner_dropped` |

## 已知局限（与论文口径保持一致）

1. LE 是**数值原型**：3 维状态是占位，未映射 LE 体系完整的 F1–F17 / E1–E7 公理；参数未标定到真实生态数据。
2. F2 的 ghost 伪力是现象学建模（副模态幅值失配 + 底噪项），尚未从界面算子不相容（提升/限制算子）解析推导。
3. F1 的吸收态性质部分来自「d 乘在通量里」这一建模选择——它是否能对应真实生殖隔离的不可逆性，需要生物学对照，**不能直接当作理论预测**。
4. 正态性/封装性判据的阈值（BC < 0.6、Var 界）目前是经验值，需要参数敏感性分析。
5. v1.6.1 的 S3/S5 结论依赖公开数据集的**当期版本**；数据集回写或修订后需重跑 `test_r3_case_upgrade.py` 复核七条预注册判据。
6. **R18 是 Python 层面的操作化，与 `formal/` 的 Lean 4 定理之间尚无对应**；
   宏库是人造的（原型 + 扰动），不是自然语言，结论是关于**机制**而不是关于语言的。
   若要进入形式化，需先写出待证命题（理论侧登记为缺口 **B-B17**）。
7. **R18 的「内约束强度」尚无有效操作化**（原 R18-b 已撤回，见 `report.md` §3.2）。
   已实现并锁进回归测试的是**二值失效的可分离性**，不是强度-反应关系。
8. `docs/xs-framework/` 的 **XS 版本号（v3.0 / v3.1）与仓库版本号（v1.6.x）是两条独立的线**：
   前者随理论主张递增，后者随代码 API / 实验 / 测试递增。文档入库**不**触发仓库主版本升级。
9. **`formal/Formal/DualConstraint.lean` 已通过 Lean 编译验证并接入**（v1.6.18，本机装 elan 后
   `lake build` 全绿、无 sorry）。原 `candidates/` 副本已按接入步骤移除。**工具链**：
   `lean-toolchain` 原指向 v4.14.0，本机下载 v4.14.0 因 SSL 连接重置反复失败，用户拍板后
   **已升级到 v4.35.0**（移除 override 复核 lake build 仍全绿），项目统一 v4.35.0。
10. **R19 / R21 是人造参数下的机制检验，不是经验验证。** R19 的四个系数
    （`e_pm` / `u_pm`）无任何经验来源；R21 的宏库是原型加扰动，不是自然语言。
    尤其：**R21 测到的全部效应都是 $O_{adm}$ 的**数量**效应**——按语义簇剪枝与
    随机剪同等数量**无差异**（原型缺乏语义结构），故该操作化应称为「$O_{adm}$ 剪枝量」，
    **不应**称为「措辞特异性」。
11. **R27 的"单峰"必须在声明的粗粒化尺度上检验。** 峰顶附近 $dU \to 0$、曲线平坦，
    未平滑时单峰检出率仅 0.045（观测噪声 $\sigma=0.15$）；移动平均窗口 $w\ge5$ 后恢复至
    $\ge0.985$。$w=21$ 是**声明式分析选择，不是拟合结果**。另：**峰位估计依赖峰后陡度**，
    B 组的 184.5 只能读作"峰位 $\ge180$"。
12. **R27 的峰位 60 / 120 / 180 是对数斜坡调度的产物，不是普适常数**；
    可迁移的是峰序关系与峰距公式的**形式**。全部结论条件于 **A_D**（B-B22 未闭合）。
13. **"单峰"在路径层是假命题。** R28 把涨落移进过程层后，路径层（$w=1,h=0$）单峰率
    为 0.000；它只在**统计层（E 层）**为真。任何对"单峰"的检验都必须在声明的
    $(w,h)$ 尺度上进行，且**优先用显著度 $h$ 而非平滑窗口 $w$**（平滑会移动峰位）。
14. **R28 的局部收敛 / 全局收敛数字（100% vs 13–34%）来自一个刻意不加停止规则的
    生灭模型**，真实过程有停止规则、会截断重开，故这些数字**不可直接迁移**；
    可迁移的是分离性本身与 $\tau_c$ 的单调性。
15. **R29 只能辨识 $\rho_{res}/\rho_{exp}$ 的比值，不能辨其绝对值。**
    绝对值的标定还须**独立知道 r 调度族**（B-B34）；目前**没有辨识调度族的程序**。
16. **定理 6.10 的峰序守恒依赖 E 与 U 共享同一 r 过程。** 若两路噪声独立，
    守恒破裂，峰序判据会以可达 **41%** 的误差率**静默判反**（B-B37）。
    且峰靠近窗口边界时会被向内拉偏，须先判定「峰是否为内点」（B-B10/T10）。
## Citation

机器可读引用见 [`CITATION.cff`](./CITATION.cff)。理论文档：《X-S-MERA 修订登记：猜想 1/2 对 A1–A10 的影响分析》v0 及后续版本；
XS 框架理论文本另见 [`docs/xs-framework/`](./docs/xs-framework/README.md)（含 v3.1 认知框架主线六篇）。

## License

代码以 **MIT License** 发布，见 [`LICENSE`](./LICENSE)。
文档中理论部分如需单独授权（例如 CC BY 4.0），请另行声明——当前仓库默认沿用 MIT。
