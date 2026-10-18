# X-S-MERA

Manifold Projection and Asymmetric Potential Well: Cognitive Dynamics Framework

> This repository contains the reference simulation code for the X‑S‑MERA theoretical framework.
> Author: [zjg]

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
| LE 单元测试 `tests/test_le_kernel.py` | ✅ `pytest tests/ -v` 全过 |
| 认知架构仿真（流形投影 + 非对称势阱，T1–T6 自证伪测试） | ❌ **尚不在本仓库**（README 此前声称的 reference code 指它，落地待补） |
| LE 公理体系 Lean4 形式化工程 `formal/` | ✅ 已通过 Lean 4.14.0 内核类型检查（wasm 运行时；native lake 见 formal/README.md） |
| LICENSE / release | ❌ 待补（徽章当前指向不存在的对象） |

## Repository structure

```
x-s-mera/
├── src/
│   ├── le_coupling_kernel.py     # LE 生命-环境耦合动力学内核（d/D、C1+C2、F1/F2、势-流分离）
│   └── kernel.py                 # （计划中）认知架构内核：快慢耦合非对称势阱
├── experiments/
│   ├── sim_le_coupling.py        # LE 耦合时序仿真：三区 / F1F2 可视化
│   ├── sim_dynamics.py           # （计划中）基础动力学
│   ├── sim_pseudo_convergence.py # （计划中）礁石伪收敛
│   ├── sim_asym_bifurcation.py   # （计划中）非对称分岔，σc 理论 -0.385 / 仿真 -0.3784
│   └── sim_scale_collapse.py     # （计划中）尺度坍塌
├── tests/
│   ├── test_le_kernel.py         # LE 内核单元测试（19 项，全过）
│   └── test_kernel.py            # （计划中）认知架构内核测试
├── docs/
│   ├── experiment_manual.md      # 实验手册（§5 LE 耦合实验）
│   └── summary_theory.md         # （计划中）理论摘要
├── formal/                       # G1-a：LE 公理体系 Lean4 形式化（L1 对象层 Horn 骨架编码）
│   ├── Formal/                   # Basic（谓词骨架）/ Axioms / Corollaries / AttackSurface
│   ├── GAPS.md                   # 编码暴露的公理缺口 B1/B2/B3/C4-full 与候选桥接
│   └── README.md                 # 构建命令与状态声明（类型检查待执行）
├── figs/
│   └── le_coupling_timeseries.png
├── requirements.txt
├── CHANGELOG.md
└── README.md
```

## Quickstart

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt        # Windows
# source .venv/bin/activate && pip install -r requirements.txt   # Linux/macOS

python experiments/sim_le_coupling.py     # 出图 + 分段统计摘要
python -m pytest tests/ -v                # 19 tests
```

## 理论 → 代码对照

| 理论条款 | 代码实现 | 检验 |
| --- | --- | --- |
| d：边界重合度标量（无量纲，非时空尺度比） | `calc_d_overlap` | `test_d_overlap_calc_range_and_monotonicity` |
| 命题 2′：C1 重合度 + C2 接口兼容，三区 | `judge_interaction_conditions` | `test_C1_C2_judge_three_region` |
| F1：无公共细化 → 直接交互禁止，通量置零 | `compute_coupling_flux` | `test_F1_branch_flux_zero` |
| F2：有重合、接口不兼容 → ghost 伪力（礁石伪收敛） | `compute_coupling_flux` F2 分支 | `test_F2_branch_ghost_flux` |
| A10：D 可封装 ⇔ 边界宏；偏离 ⇔ M-不完备 | `LEProjectionSystem.test_D_encapsulation` | `test_D_encapsulation_*` |
| A11-R：势只调制增益，不进入状态方程 | 通量公式中势仅出现在 gain | `test_potential_never_enters_state_equation` |

## 已知局限（与论文口径保持一致）

1. LE 是**数值原型**：3 维状态是占位，未映射 LE 体系完整的 F1–F17 / E1–E7 公理；参数未标定到真实生态数据。
2. F2 的 ghost 伪力是现象学建模（副模态幅值失配 + 底噪项），尚未从界面算子不相容（提升/限制算子）解析推导。
3. F1 的吸收态性质部分来自「d 乘在通量里」这一建模选择——它是否能对应真实生殖隔离的不可逆性，需要生物学对照，**不能直接当作理论预测**。
4. 正态性/封装性判据的阈值（BC < 0.6、Var 界）目前是经验值，需要参数敏感性分析。

## Citation

见 `CITATION.cff`（待补）。理论文档：《X-S-MERA 修订登记：猜想 1/2 对 A1–A10 的影响分析》v0 及后续版本。
