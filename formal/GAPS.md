# GAPS：Γ_LE 推论缺口与候选桥接公理（G1-a 的直接产出）

G1-a 的形式化价值不止于"无矛盾"——编码过程把归档版第 7 节
"推论全量可演绎导出，非公理"的声明变成了可核验命题，并暴露了缺口。
以下推论在 Horn skeleton 层**不可导出**；候选桥接公理入 Γ 前
须经修订登记走正式修订周期，并**重跑 AttackSurface 的一致性检查**。

## B1（归档推论 1：衰老）

- 归档文本：结构持续损伤累积 + 能量修复上限约束 ⇒ 必然产生个体老化。
- 缺口：F4 只断言"损伤不可逆累积"；F9 只断言"能量中断加速损伤"。
  无任何公理断言"损伤必然突破阈值"，也无"修复存在上限"的谓词。
- 候选桥接（任选其一，须论证）：
  - B1a：引入 `EventuallyDamageOverThresh : Sys → Prop` + 公理
    `LifeS s → DamageAccIrrev s → EventuallyDamageOverThresh s`
  - B1b：引入 `RepairCap : Sys → Prop` + 公理
    `DamageAccIrrev s → RepairCap s → EventuallyDamageOverThresh s`
- 若两者都不采纳：衰老推论降级为"机制存在性"
  （`AgingMechanism := DamageAccIrrev`，F4 直接给出）而非"必然性"。
  该降级须登记——与 R5 的"逐语境宏 M"降级机制同构。

## B2（归档推论 3：演化）

- 归档文本：无向变异提供素材 + 资源环境约束筛选 ⇒ 定向演化。
- 缺口：F14 给变异、F13 给竞争，无公理连接二者到"定向演化"
  （Γ_LE 中甚至没有 `DirectedEvolution` 谓词）。
- 候选桥接：引入 `DirectedEvolution : Sys → Prop` + 公理
  `UndirectedHeritableVariation s → Competition s → DirectedEvolution s`。
- **前置条件**：R5 已证"定向演化"是逐语境命题
  （Ns≳5 语境成立 / Ns→0 不成立）。B2 若入 Γ 必须带 C_q 资格条件
  （尺度/参数区间申报，R5 字段），否则与 R5 结论直接冲突。

## B3（归档推论 6：吸引子分类）

- 归档文本：体系自发抵达四类终态（稳态 / 周期振荡 / 生态锁定 / 全局灭绝）。
- 缺口：终态分类是动力学稳定性陈述，需要状态空间与流结构，
  FOL 谓词层（Horn skeleton）结构性不可达。
- 处置：**不在 Γ_LE 增补公理**；登记为 L0 层目标，
  交仿真层核验（x-s-mera 内核的四类终态相图测绘，属"动态仿真"遗留事项）。

## C4-full（归档推论 4：耗散结构，完整版）

- 骨架版可证（Corollaries.lean：F6→F8、F9 前件）。
- 完整版"生命以环境熵增为代价维持有序"须 M1 + 一个未在签名中的前提：
  生命-环境**合系统**的孤立性（Isolated 作用于合系统）。
- 候选处置：扩展签名为 `Isolated : Sys → Sys → Prop`（两元合系统孤立性），
  或登记完整版为 M1 的 L2 操作层实例（仿真验证熵账本闭合），不进 L1。

## 相对一致性注记

Γ_LE 当前片段（正 Horn + coupling_iff + 等价关系性）的一致性
由元论证支持：无否定公理 ⇒ 无矛盾源。该论证在引入 B1/B2
（含否定或新谓词）后可能失效，届时须改用模型构造路线
（在 Lean 中对 Sys 做平凡解释并验证公理全真——
注意该模型同时证明"仅平凡模型存在"的风险，即 G2 的非平凡性
必须靠 L2 层见证（LE 内核）而非 Lean 模型——两层分工见 v3.3.1）。


## C1（认知双约束：XS v3.1 补篇三 / R18 的形式化）

- **状态：已接入（2026-10-24，CHANGELOG v1.6.18）。** 规格现位于
  [`Formal/DualConstraint.lean`](./Formal/DualConstraint.lean)，已加入 `Formal.lean` 的 import，
  `lake build` 编译通过、无 sorry。原 `candidates/DualConstraint.lean` 已按"移入"移除。
- **接入记录（诚实）**：`lean-toolchain` 原指向 `leanprover/lean4:v4.14.0`，接入机从
  `releases.lean-lang.org` 下载 v4.14.0 反复因 **SSL 连接被重置** 失败（stable v4.35.0 可下载）；
  先用 `elan override` 以 v4.35.0 编译验证通过（Lean 4 基础语法兼容），**随后用户拍板将
  `lean-toolchain` 升级到 v4.35.0**，移除 override 复核 `lake build` 仍全绿，项目统一 v4.35.0。
- **接入步骤（已完成）**：单独编译 → 移入 `Formal/DualConstraint.lean` →
  `Formal.lean` 加 `import Formal.DualConstraint` → CHANGELOG 登记 → 本条标记"已接入"。

### 形式化范围（受 XS 侧未闭合缺口的限制）

- **已编码**：两类失效的**可分离性**（DC3a/DC3b）与**判别签名 R18-a**
  （DC4 外失效 ⇒ σ 不降；DC5 内失效 ⇒ σ 降），以及 T1–T6 六条可导出定理。
  这些都在**操作层**，不依赖 XS 侧尚未闭合的 B-B19。
- **结构性不可达（登记为 X1，与上面 B3 同构）**：R18 的**定量签名**
  —— `|σ(外失效) − σ(base)| ≤ 0.006` 而 `σ(base) − σ(内失效) ≥ 0.199`，
  且二者相差一个数量级。这是数值不等式，Horn skeleton 不含算术结构，**FOL 层不可表达**。
  处置：不在本层增补公理，由 `experiments/r18_dual_constraint/` 的两个回归测试锁死。
- **暂不编码（登记为 X2）**：XS 侧 B-B19（动力学层 Ψ 与信息层 U 的桥）未闭合，
  因此「收敛 ⟺ U ≈ 0」「定理 5.2′」「推论 5.3（酝酿 / 反刍判别）」**暂不进入形式化**。
  它们需要 `Ctx → Rat` 之类的数值映射，会**改变本层签名**，须单独立项
  （建议 `Formal/CognitiveBridge.lean`）。

### 公理地位（须走修订周期）

DC3–DC8 是**经验结论的公理化**，不是逻辑真理。其中 DC6 / DC7（两类失效各自主导一种
失败模式）来自 800 次/臂的实测。若出现反例，须降级或加 C_q 资格条件——
与上面 B2 要求资格条件的处置同构。

### 跨层耦合的前置条件

本文件的认知层与 Γ_LE 层**不共享任何谓词**，因此不构成对 Γ_LE 的扩展。
若日后要把「理解者」建模为 `Sys`（从而认知层与生命-环境层耦合），
**必须**重跑 `AttackSurface.lean` 的一致性检查——那是一次真正的跨层合并，
不能只靠加 import。
