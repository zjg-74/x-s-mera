/-
# DualConstraint：认知双约束（XS v3.1 补篇三 / R18）的 Horn skeleton 编码 —— **候选，未编译**

## 状态声明（务必先读）

本文件**未经 Lean 4 编译验证**。生成它的机器上没有 Lean 工具链
（`lean-toolchain` 指定 `leanprover/lean4:v4.14.0`，但 `lean` / `lake` 均不可用）。
因此它被放在 `formal/candidates/` —— 该目录不在 `lakefile.toml` 的 `Formal`
目标下，**不会被 `lake build` 编译**，不会破坏既有构建。

**接入步骤（在有工具链的机器上执行）**：
1. `lean Formal/candidates/DualConstraint.lean` 单独编译验证（或临时移入 `Formal/` 后 `lake build`）；
2. 语法 / 证明修正后，移入 `formal/Formal/DualConstraint.lean`；
3. 在 `Formal.lean` 增加 `import Formal.DualConstraint`；
4. 在 `CHANGELOG.md` 登记，并把 `GAPS.md` 中 C1 条目标记为"已接入"。

## 形式化的范围（重要限制）

XS 理论侧的缺口 **B-B19**（信息层 `U` 与动力学层 `Psi` 的桥）**尚未闭合**。
因此本文件**只编码不依赖该桥的部分**：两类失效的**可分离性**与
**判别签名 R18-a**。二者都在操作层，不涉及 `U` / `Psi`。

**结构性不可达的部分**（与 `GAPS.md` 中 B3「吸引子分类」同构，交仿真层）：
R18 的**定量签名** —— `|sigma(外失效) - sigma(base)| <= 0.006` 而
`sigma(base) - sigma(内失效) >= 0.199`，且二者相差一个数量级。
这是**数值**断言，Horn skeleton 的 FOL 层无法表达。

## 公理地位声明

DC3–DC7 是**经验结论的公理化**，不是从更基本命题推出来的。
本仓库 `Axioms.lean` 中的 F1–F16 本就是「事实公理」，做法一致；
但按修订纪律，**新增公理须走正式修订周期**，并注意：

- DC6/DC7（两类失效各自主导一种失败模式）来自
  `experiments/r18_dual_constraint/` 的 800 次/臂实测，**不是逻辑真理**；
- 若后续出现反例（某系统外失效却 F1 主导），DC6 须降级或加资格条件，
  与 `GAPS.md` 中 B2 要求 C_q 资格条件的处置同构。

## 与既有文件的关系

本文件不引入新的论域类型，复用 `Formal.Basic` 的风格（谓词 + Horn 子句）。
认知层与 LE（生命–环境）层**不共享谓词**，本文件不构成对 Γ_LE 的扩展，
因此**不触发** `AttackSurface.lean` 的一致性检查范围之外的风险。
若日后要把认知层与 LE 层耦合（例如把「理解者」建模为 Sys），
**必须**重跑 AttackSurface 检查 —— 那是一次真正的跨层合并。
-/

namespace XSMERA

-- ========== 论域 ==========

/-- 认知情境：一次「理解这句话」的完整设置（输入 + 任务 + 理解者） -/
axiom Ctx : Type

-- ========== 约束的在位性 ==========

/-- 外约束在位：输入的物理形态未被替换（对象未被换掉） -/
axiom OuterActive : Ctx → Prop

/-- 内约束在位：任务仍是「理解」，未被替换为续写 / 改写 / 批评等操作 -/
axiom InnerActive : Ctx → Prop

-- ========== 失败模式（R18 的操作化定义） ==========

/-- F1：未加载 / 不收敛（无可用宏，或达不到接受条件） -/
axiom F1dominant : Ctx → Prop

/-- F2：伪收敛（加载成功且满足落定条件，但事后验证为错） -/
axiom F2dominant : Ctx → Prop

-- ========== 判别签名 R18-a：终局加载度 sigma 的行为 ==========

/-- 终局 sigma 相对基线**未下降**（外约束是方向约束：错了但一样确信） -/
axiom SigmaUndropped : Ctx → Prop

/-- 终局 sigma 相对基线**显著下降**（内约束是强度约束：不够确信就停） -/
axiom SigmaDropped : Ctx → Prop

-- ========== 内约束强度（B-B21 的操作化：O_adm 剪枝量） ==========

/-- 内约束对允许操作集做了剪枝 -/
axiom AdmPruned : Ctx → Prop

-- ========== 结构性公理 ==========

/-- 加载函数（sigma 的计算方式）未被改动 —— DC1 的后件 -/
axiom LoadFnPreserved : Ctx → Prop

/-- 停止规则被改动 —— DC2 的后件 -/
axiom StopRuleChanged : Ctx → Prop

/-- DC1 外约束作用于**对象**：失效只替换输入，不改变加载函数与停止规则。 -/
axiom DC1 : ∀ c : Ctx, ¬ OuterActive c → LoadFnPreserved c

/-- DC2 内约束作用于**操作集与停止规则**：失效改变任务，不改变输入。 -/
axiom DC2 : ∀ c : Ctx, ¬ InnerActive c → StopRuleChanged c

/-- DC3a 逻辑独立性之一：存在「外失效而内未失效」的情境。 -/
axiom DC3a : ∃ c : Ctx, ¬ OuterActive c ∧ InnerActive c

/-- DC3b 逻辑独立性之二：存在「内失效而外未失效」的情境。 -/
axiom DC3b : ∃ c : Ctx, OuterActive c ∧ ¬ InnerActive c

/-- DC4 外失效 ⇒ 终局 sigma 不下降（R18-a 的上半支）。 -/
axiom DC4 : ∀ c : Ctx, ¬ OuterActive c → SigmaUndropped c

/-- DC5 内失效 ⇒ 终局 sigma 显著下降（R18-a 的下半支）。 -/
axiom DC5 : ∀ c : Ctx, ¬ InnerActive c → SigmaDropped c

/-- DC6 外失效 ⇒ F2 主导（对象被替换 ⇒ 加载了不来自这句话的东西）。 -/
axiom DC6 : ∀ c : Ctx, ¬ OuterActive c → F2dominant c

/-- DC7 内失效 ⇒ F1 主导（任务被替换 ⇒ 不收敛）。 -/
axiom DC7 : ∀ c : Ctx, ¬ InnerActive c → F1dominant c

-- ========== 互斥性（结构性，非经验） ==========

/-- 同一情境的 sigma 不可能同时「未下降」又「显著下降」。 -/
axiom SigmaExcl : ∀ c : Ctx, SigmaUndropped c → ¬ SigmaDropped c

/-- 同一情境不可能同时 F1 主导与 F2 主导。 -/
axiom F1F2Excl : ∀ c : Ctx, F1dominant c → ¬ F2dominant c

/-- 内约束在位 ⇒ 允许操作集被剪枝（B-B21：Q = log2(|O_all|/|O_adm|) > 0）。 -/
axiom DC8 : ∀ c : Ctx, InnerActive c → AdmPruned c

-- ========== 可导出定理 ==========

/-- T1 外失效不导致 F1 上升：外失效给出 F2，而 F2 与 F1 互斥。 -/
theorem T1_outer_failure_not_F1 (c : Ctx) (h : ¬ OuterActive c) : ¬ F1dominant c := by
  intro hf1
  exact (F1F2Excl c hf1) (DC6 c h)

/-- T2 内失效导致 F1 主导（DC7 的直接重述）。 -/
theorem T2_inner_failure_is_F1 (c : Ctx) (h : ¬ InnerActive c) : F1dominant c := by
  exact DC7 c h

/-- T3 判别签名 R18-a：两类失效的 sigma 行为**可分离**。
    注意这里证明的是**可分离性**（两个签名不同），不是数值差距 ——
    数值差距（相差一个数量级）在 FOL 层不可表达，见文件头的状态声明。 -/
theorem T3_sigma_signature_separable (c1 c2 : Ctx)
    (h1 : ¬ OuterActive c1) (h2 : ¬ InnerActive c2) :
    SigmaUndropped c1 ∧ SigmaDropped c2 := by
  exact ⟨DC4 c1 h1, DC5 c2 h2⟩

/-- T4 单侧失效不会被另一侧"吞掉"：外失效且内约束在位时，仍然是 F2 主导。
    （对应 R18 的 A 臂：只篡改输入，任务仍为「理解」。） -/
theorem T4_outer_only_is_F2 (c : Ctx)
    (ho : ¬ OuterActive c) (_hi : InnerActive c) : F2dominant c := by
  exact DC6 c ho

/-- T5 内失效且外约束在位时，仍然是 F1 主导。
    （对应 R18 的 B 臂：输入正常，任务改为「续写」。） -/
theorem T5_inner_only_is_F1 (c : Ctx)
    (_ho : OuterActive c) (hi : ¬ InnerActive c) : F1dominant c := by
  exact DC7 c hi

/-- T6 双失效时 sigma 仍下降 —— 内约束主导 sigma。
    （对应 R18 的 C 臂：sigma 与 B 臂同为 0.375，与 base 的 0.646 显著不同。）
    这是 R18 的一个**非平凡**观测：内约束对 sigma 的作用优先于外约束。 -/
theorem T6_inner_dominates_sigma (c : Ctx)
    (_ho : ¬ OuterActive c) (hi : ¬ InnerActive c) : SigmaDropped c := by
  exact DC5 c hi

-- ========== 不可达登记（不编码，仅记录） ==========

/-
## X1（FOL 层不可达，交仿真层）

R18 的定量签名：

    |sigma(外失效) - sigma(base)| <= 0.006
    sigma(base) - sigma(内失效)   >= 0.199
    且后者 / 前者 >= 5（相差一个数量级）

这是**带数值的**不等式，本文件的 Horn skeleton 不含算术结构，
结构性不可表达。处置与 `GAPS.md` 中 B3（吸引子分类）同构：
**不在本层增补公理**，由 `experiments/r18_dual_constraint/` 的
回归测试 `test_signature_sigma_outer_unchanged_inner_dropped` 与
`test_signature_holds_at_higher_resolution` 锁死。

## X2（依赖未闭合的 B-B19，暂不编码）

XS 理论侧的「动力学层 Psi」与「信息层 U」的桥尚未闭合
（见 `docs/xs-framework/06-*.md` 补篇五 §1）。因此以下命题**暂不进入形式化**：

  * U 的单调性与三态的对应（收敛 ⟺ U ≈ 0）
  * 定理 5.2'（sign(psi_bar) = -sign(dU/dk) ⟺ rho_res > rho_exp）
  * 推论 5.3（震荡期 U 的走向判别酝酿 / 反刍）

B-B19 闭合后，这些须单独立项（建议文件 `Formal/CognitiveBridge.lean`），
因为它们需要 `Ctx → Rat` 之类的数值映射，会**改变本层的签名**。
-/

end XSMERA
