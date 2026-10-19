# XS 框架理论文档集

本目录收录 **XS 框架**（尺度坍缩与多投影统一）的理论文本，与仓库 `formal/`（Lean 4）
和 `src/`（Python 数值原型）构成三层：**理论层 / 形式层 / 数值层**。

版本线说明（重要）：

- **XS 框架版本号（v3.0 / v3.1）与仓库代码版本号（v1.6.x）是两条独立的线**。
  仓库版本遵循语义化版本，仅在代码 API / 实验 / 测试变动时递增；
  XS 版本在理论主张变动时递增。**本文集入库本身不触发仓库主版本升级。**
- 文档按 **U1（取代不删除）** 纪律维护：后出的补篇**不修改**前篇正文，
  只在前篇顶部插入条件化标注互指。因此**出现"同一命题有两个版本"是预期的**，
  以**编号最大者**为准，前篇保留作为退化情形的兜底与审计线索。

## 一、v3.1 认知框架主线（2026-10-19，本轮新增）

从**认知二分**出发，经**宏加载**得到 XSM，讨论认知如何**三分（收敛 / 发散 / 震荡）**，
并把认知过程类比为**基于事实公理与运算公理的系统演绎**。

| 文件 | 内容 | 关键装置 |
| --- | --- | --- |
| [`01-v3.1-main-dichotomy-macro-loading-trichotomy.md`](./01-v3.1-main-dichotomy-macro-loading-trichotomy.md) | **主文**。§2 认知二分（X = 象 / S = 运算）→ §3 事实公理 F 与运算公理 O → §4 宏加载 → §5 三分 → §6 概率三分 → §7 与既有体系对齐 → §8 缺口 → §9 可检验推论 R1–R4 | 加载度 σ；编译型 / 调用型加载；**三条偏离 → 三个归宿**（事实开放→概率性；资源有界→发散；承诺可撤销→伪收敛） |
| [`02-v3.1-supp1-chaos-boundary-negation-closure.md`](./02-v3.1-supp1-chaos-boundary-negation-closure.md) | **补篇一**：混沌边界与否定闭合 | 定理 3.2（非单调性局部化到 `eff`）；极性定理（**混合极性 = 震荡的动力学来源**）；σ 改为基于 `eff` 的可废后果 |
| [`03-v3.1-supp2-potential-drift-diffusion.md`](./03-v3.1-supp2-potential-drift-diffusion.md) | **补篇二**：势的漂移-扩散重构 | **势 Ψ := P − Q**（双约束下两股反向力之净差）；**Ψ = ψ̄ + ξ**（统计定向 + 局部回落）；定理 6.2″（ψ̄ ≠ 0 ⇒ 几乎必然触及边界） |
| [`04-v3.1-supp3-dual-constraint-PQ.md`](./04-v3.1-supp3-dual-constraint-PQ.md) | **补篇三**：双约束的构成与 P/Q 构造 | 外约束 = 输入的物理约束（构成性，对 X）；内约束 = 必须"理解"（规范性，对 S）；**两类失效分工**（外失效→F2/幻觉，内失效→F1/不收敛） |
| [`05-v3.1-supp4-calibration-commensurability.md`](./05-v3.1-supp4-calibration-commensurability.md) | **补篇四**：标定、通约性自我修正、B-B14 裁决 |
| [`06-v3.1-supp5-bridge-constraint-strength.md`](./06-v3.1-supp5-bridge-constraint-strength.md) | **补篇五**：U 与 Ψ 的桥（B-B19）、内约束强度的操作化（B-B21）、P 的理解者依赖（B-B16）、R18 的 Lean 命题（B-B17） | **定理 5.2′**（猜想 5.2 无条件形式**被证否**，带 $ho_{	ext{res}} > ho_{	ext{exp}}$ 条件成立）；**推论 5.3**（震荡期 U 的走向判别**酝酿 vs 反刍**）；内约束强度 = $O_{	ext{adm}}$ 剪枝量；**主文 §5.4 成因 2 修正**：过早锁定的危害是"正确选项不在 $A_k$ 中"，不是"候选太少" | **Q := I(M; C_in \| φ)**（取代不可通约的原定义）；**U := P − Q = H(M \| φ, C_in)**（净未决量，非负）；P 的标定程序（经验分布熵 + Miller–Madow 修正） |

### 三分的成立条件（当前辖域，逐篇收紧）

| 版本 | 条件 |
| --- | --- |
| 主文 §1.2 | （无条件主张） |
| 补篇一 §7 | 资源有界 `B < ∞` ∧ 未决代价 `μ > 0` |
| **补篇二 §8（现行）** | **资源有界 `B < ∞` ∧ 存在统计定向 `ψ̄ ≠ 0`** |

`ψ̄ = 0` 时熵轨迹为鞅，震荡与混沌可持久 → **四分**，且吸引盆分形化使概率权重不可鲁棒赋值。
**局部回落（`Var(ξ) > 0`）不妨碍归宿，只延迟归宿。**

### 认知存在的条件（补篇三/四，另一层，勿与上表混淆）

`0 < H(M|φ) < ∞` ∧ `Q > 0`。`P = 0` → 理解退化为解码；外约束消失 → 理解退化为生成（幻觉）。

## 二、v3.0 统一版与元规则

| 文件 | 内容 |
| --- | --- |
| [`10-v3.0-unified-framework-paper.md`](./10-v3.0-unified-framework-paper.md) | 统一重构版论文：§2 方法论（基底-完备七门、S1/S2/S3、P0–P8）→ §3 XS 体系 → §4–6 三个应用（量子-经典、MMS、观察者诊断）→ §7 工程反哺 → §9 U1–U4 |
| [`11-three-version-unification-and-engineering-feedback.md`](./11-three-version-unification-and-engineering-feedback.md) | 三版（8 月 / 10 月 / 工程）术语台账 + 四条反哺通道 A–D + 推论 R/R′/R″/R‴ |
| [`12-meta-rules-update-discipline-and-corollary-P.md`](./12-meta-rules-update-discipline-and-corollary-P.md) | **U1–U4 更新纪律**（U1 取代不删除 / U2 类型 E·I·X+·X− / U3 依赖追踪 / U4 冻结窗口）+ Lakatos 进步判据 + 推论 P |
| [`13-decisions-manifold-adjudication.md`](./13-decisions-manifold-adjudication.md) | 三项决定：流形归属**裁决为结构层**；替换 D；三案例重论 |

## 三、8 月原始稿与审计（存档）

| 文件 | 内容 |
| --- | --- |
| [`20-august-six-articles-deduped.md`](./20-august-six-articles-deduped.md) | 8 月六篇正文去重版（原文件三重重复 + 段内重复，已清理；流形替换 A–I 已应用） |
| [`21-audit-august-return-to-mainline.md`](./21-audit-august-return-to-mainline.md) | 回归主线审计：①方法论 → ②体系 → ③应用 |
| [`22-audit-methodology-first.md`](./22-audit-methodology-first.md) | 方法论优先路径 |
| [`23-audit-basis-vector-space-metaphor.md`](./23-audit-basis-vector-space-metaphor.md) | **"n 维向量空间 + n 个基底向量即事件公理"隐喻的审计**：演绎闭包不满足 Steinitz 交换性（反例：`p ∈ Cn({p∧q})`、`p ∉ Cn(∅)`，但 `p∧q ∉ Cn({p})`）→ Cn 不是拟阵 → 无维数不变性；正确量是 `rank(π)` |
| [`24-audit-linear-nonlinear-scale-evolution.md`](./24-audit-linear-nonlinear-scale-evolution.md) | "X1 线性 → 高层非线性"命题审计：三个反例（GR / Yang–Mills / NS 在底层即非线性；Onsager / Kubo 在宏观仍线性；Ising 双线性但 Lee–Yang 奇点） |
| [`25-audit-nonprofessional-observer.md`](./25-audit-nonprofessional-observer.md) | 非专业观察者言论重估：Hesse 中性类比 / Gentner 结构映射 / Chi–Feltovich–Glaser / Bachelard 认识论障碍 |

## 四、缺口登记（跨文件，当前未闭合项）

按优先级排列。编号 `B-Bn` 为 v3.1 认知线专用，与 8 月/10 月其他缺口表分开计数。

| 编号 | 优先级 | 内容 | 出处 |
| --- | --- | --- | --- |
| **B-B19** | **中（已定位，未闭合）** | 猜想 5.2 作为**无条件命题已被证否**（三个反例模型并存）；在 A_D 下成为**定理 5.2′**，条件 $ho_{	ext{res}} > ho_{	ext{exp}}$（可测，见 R23）。**跨层引用由无条件禁令改为条件化许可** | 补篇五 §1 |
| **B-B22** | **高** | **A_D（展开/消解双速率线性分解）本身是未证的建模假设**。定理 5.2′ 与推论 5.3 全部依赖它。证否签名：**同一理解者同时出现 $(+,-)$ 与 $(-,+)$**（R24 未跑） | 补篇五 §1.2、§5 |
| **B-B11** | 高 | 主文 §3.5-B 缺"资源耗尽 → `A_k ≠ ∅` → 伪收敛"这一支 | 补篇二 §3 |
| **B-B17** | 中（部分推进） | 命题规格已写出：`formal/candidates/DualConstraint.lean`（DC1–DC8 + T1–T6 + 不可达登记 X1/X2）。**本机无 Lean 工具链，未编译验证**，故置于非编译目录；见 `formal/GAPS.md` 条目 C1 | 补篇五 §4 |
| **B-B16** | 已闭合 | 已正式改写为 $P := H(\mathcal{M}\midarphi,\ \Omega)$（$\Omega$ = 理解者）。**不仅 U 依赖理解者，连"两层能否互引"也依赖理解者**（推论 3.1）。**新增 B-B25** | 补篇五 §3 |
| **B-B23** | 中 | $ho_{	ext{res}}$ 与 $ho_{	ext{exp}}$ 的**单独标定**未给（推论 5.3 只给二者**大小关系**） | 补篇五 §5 |
| **B-B24** | 中 | 原型缺乏语义结构：按簇剪枝 ≈ 随机剪枝，**无法检验"措辞"的结构效应**。故该操作化应称「$O_{	ext{adm}}$ 剪枝量」，**不应**称「措辞特异性」 | 补篇五 §2.6 |
| **B-B25** | 中 | $\hat P$ 的跨理解者混合估计**系统性低估** $P$，偏倚 $= I(\mathcal{M};\Omega\midarphi)$，未给修正程序。此前 R15′/R20 的数值**不可跨群体比较** | 补篇五 §3.3 |
| **B-B12** | 低（部分闭合） | P 的数值标定：程序已给（经验分布熵 + Miller–Madow），实证估计待做 | 补篇四 §3 |
| **B-B14** | 已裁决 | 三组二分（X/S、F/O、外/内）**只按弱形式保留**，不主张严格同构 | 补篇四 §1 |
| **B-B21** | **已闭合** | 内约束强度 = $O_{	ext{adm}}$ 剪枝量 $Q = \log_2(\lvert O_{	ext{all}}vert/\lvert O_{	ext{adm}}vert)$（不依赖 $	au$ 或注入策略）。R21 实测 F2 呈**倒 U**（峰值 $Qpprox1.0$）。**注意**：本人提出的"虚假唯一性"预言**被证否**，由此修正主文 §5.4 成因 2 | 补篇五 §2 |

已闭合：B-B6（三个"势"命名消歧）、B-B13（Q 的通约性）、B-B10（双约束构成，主要部分）。

## 五、与代码层的对应

| 理论条款 | 代码 / 形式化 | 状态 |
| --- | --- | --- |
| 两条失效分工（外→F2、内→F1） | `experiments/r18_dual_constraint/` | **已跑通**，见该目录 `report.md` |
| A11-R 语义势只调增益、不进状态方程 | `test_potential_never_enters_state_equation` | 既有（134 项之一） |
| 加载度 σ ↔ 工程版 d（边界重合度） | — | 登记为推论 R1，**尚未实现** |
| F1 ↔ 发散、F2 ↔ 伪收敛 | — | 登记为推论 R2 |
| 双约束的可分离失效 | — | 登记为推论 R18 |

**R18 的判决性实验建议（下一步）**：R8（扫语义势增益看是否出现 logistic 式分岔级联，
检验"势改步长"的合法作用）与 R11（人为违反 A11-R，看是否产出高加载度但事后为错，
检验"势越权"的后果）可在现有离散迭代原型上直接执行，与 R18 构成一组三条。

## 六、读法建议

1. 首次接触：主文 `01` → 补篇二 `03`（势）→ 补篇三 `04`（双约束）。
2. 只关心可检验结论：主文 §9（R1–R4）+ 补篇一 §6（R5–R8）+ 补篇二 §6（R9–R12）
   + 补篇三 §6（R15–R18）+ 补篇四 §6（R19–R20）+ 补篇五 §6（R23–R26）。
3. 审计/复核：先看本文件第四节缺口表，再看各篇的"诚实收缩"小节。
