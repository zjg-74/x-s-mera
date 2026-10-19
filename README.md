# X-S-MERA

Manifold Projection and Asymmetric Potential Well: Cognitive Dynamics Framework

> This repository contains the reference simulation code for the X‑S‑MERA theoretical framework.
> Author: [zjg]
>
> Current version: **v1.6.3**（2026-10-19） · Tests: **155 passed** · License: **MIT**

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
| 认知架构仿真（流形投影 + 非对称势阱，T1–T6 自证伪测试） | ❌ **尚不在本仓库**（README 此前声称的 reference code 指它，落地待补） |
| LE 公理体系 Lean4 形式化工程 `formal/` | ✅ 已通过 Lean 4.14.0 内核类型检查（wasm 运行时；native lake 见 formal/README.md） |
| LICENSE | ✅ MIT（见 `LICENSE`）；release 待发布 |
| CITATION | ✅ `CITATION.cff` 已补 |
| **XS 框架理论文档集** `docs/xs-framework/`（v3.1 认知框架主线 + v3.0 统一版 + 8 月稿与审计） | ✅ 已入库（v1.6.3），**16 篇**；索引见该目录 `README.md` |
| **R18 双约束失效可分离性实验** `experiments/r18_dual_constraint/` | ✅ 已实现（v1.6.2），7 项测试；报告见该目录 `report.md` |
| **R19 U/Psi 桥一致性检验** `experiments/r19_bridge/` | ✅ 已实现（v1.6.3），7 项测试 |
| **R21 内约束强度 = O_adm 剪枝量** `experiments/r21_constraint_specificity/` | ✅ 已实现（v1.6.3），7 项测试 |
| **R18 的 Lean 命题规格** `formal/candidates/DualConstraint.lean` | ⚠️ **候选，未编译验证**（本机无 Lean 工具链），置于非编译目录；见 `formal/GAPS.md` 条目 C1 |

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
│   ├── r18_dual_constraint/      # R18：双约束失效可分离性（5 脚本 + 5 json + report.md）
│   ├── r19_bridge/               # R19：U/Psi 桥的双速率分解一致性检验（1 脚本 + json + report.md）
│   └── r21_constraint_specificity/  # R21：内约束强度 = O_adm 剪枝量（3 脚本 + 3 json + report.md）
├── tests/                        # 155 项单元测试，全过（14 个测试文件）
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
python -m pytest tests/ -v                # 155 tests
python experiments/r18_dual_constraint/r18_dual_constraint.py   # R18 主实验（仅标准库）
python experiments/r19_bridge/r19_bridge.py                     # R19（仅标准库）
python experiments/r21_constraint_specificity/r21_constraint_specificity.py  # R21（仅标准库）
```

复现须知：

- 需要 Python ≥ 3.10（依赖 `numpy>=2.0`）。
- `tests/test_r3_case_upgrade.py` 含数据完整性校验（防回写漂移），耗时明显高于其余文件；
  只做快速自检可先跑 `python -m pytest tests/ -v -x --ignore=tests/test_r3_case_upgrade.py`。
- 全套 141 项的分文件分布：le_kernel 19 / a10 16 / a9 13 / r3_case_upgrade 14 / r4 12 /
  t1d 12 / r3_sufficiency 11 / relay_chain 10 / tp 10 / a0_bootstrap 9 / d_form_audit 8 /
  **r18_dual_constraint 7**（v1.6.2）/ **r19_bridge 7** + **r21_constraint_specificity 7**（v1.6.3）。
- `experiments/r18_dual_constraint/` 下的脚本**只依赖标准库**，不需要 numpy；
  运行时把 json 写进当前工作目录，故请在该目录内执行（见 `report.md` §7）。

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
9. **`formal/candidates/DualConstraint.lean` 未经 Lean 编译验证**（本机无 `lean` / `lake`）。
   它故意放在 `candidates/`（不在 `lakefile.toml` 的 `Formal` 目标下，不会被编译），
   以免破坏构建。验证后移入 `Formal/` 并 import 即完成接入，步骤见 `formal/GAPS.md` 条目 C1。
10. **R19 / R21 是人造参数下的机制检验，不是经验验证。** R19 的四个系数
    （`e_pm` / `u_pm`）无任何经验来源；R21 的宏库是原型加扰动，不是自然语言。
    尤其：**R21 测到的全部效应都是 $O_{adm}$ 的**数量**效应**——按语义簇剪枝与
    随机剪同等数量**无差异**（原型缺乏语义结构），故该操作化应称为「$O_{adm}$ 剪枝量」，
    **不应**称为「措辞特异性」。

## Citation

机器可读引用见 [`CITATION.cff`](./CITATION.cff)。理论文档：《X-S-MERA 修订登记：猜想 1/2 对 A1–A10 的影响分析》v0 及后续版本；
XS 框架理论文本另见 [`docs/xs-framework/`](./docs/xs-framework/README.md)（含 v3.1 认知框架主线六篇）。

## License

代码以 **MIT License** 发布，见 [`LICENSE`](./LICENSE)。
文档中理论部分如需单独授权（例如 CC BY 4.0），请另行声明——当前仓库默认沿用 MIT。
