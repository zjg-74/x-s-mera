/-!
# Corollaries：归档版第 7 节六条推论的演绎审计

归档版声称"推论全量可演绎导出，非公理"。本文件把该声明变成可核验命题：
能证的给出证明；证不出的不虚构证明，登记到 GAPS.md（桥接公理候选）。

**本文件纪律：禁止使用 sorry/admit。** 证不出的推论以注释登记，不用占位符伪造。
-/
import Formal.Basic
import Formal.Axioms

namespace XSMERA

-- ========== 第 2 节后件（补充审计项，非归档第 7 节编号） ==========

/-- 尺度不一致 ⇒ 不能双向稳定耦合。由核心定理 iff 的逆否直接导出。 -/
theorem mismatch_forbids_coupling :
    ∀ l e : Sys, LifeS l → EnvS e → ¬ SameScale l e → ¬ Coupled l e := by
  intro l e _ _ hns hc
  exact hns ((coupling_iff l e).mp hc)

-- ========== 推论 1（衰老）——缺口登记 ==========
-- 归档文本：损伤持续累积 + 能量修复上限约束 ⇒ 必然产生个体老化。
-- 审计结论：F4 只断言"累积不可逆"，F9 只断言"能量中断加速损伤"；
--   无任何公理断言"损伤必然突破阈值"或"修复存在上限"。
--   "必然老化"在 Horn skeleton 层不可导出 → 候选桥接公理 B1（见 GAPS.md）。
-- 可证的部分只是"老化机制存在"（= F4 重述）：
theorem aging_mechanism_exists :
    ∀ s : Sys, LifeS s → DamageAccIrrev s := fun s hs => F4 s hs

-- ========== 推论 2（双层消亡）——公理重述登记 ==========

/-- 个体层：损伤破阈值 ⇒ 死亡。即 F5 本身。
    审计登记：归档称"推论"，但个体层无独立演绎内容——它是公理重述。 -/
theorem individual_death :
    ∀ s : Sys, LifeS s → DamageOverThresh s → Death s := fun s hs hd => F5 s hs hd

/-- 谱系层：增殖率低于替换阈值 ⇒ 谱系萎缩至灭绝。即 F15 本身（重述登记同上）。 -/
theorem lineage_extinction :
    ∀ s : Sys, LineageS s → BelowReplacementRate s → LineageDeclineToExtinction s :=
  fun s hs hb => F15 s hs hb

-- ========== 推论 3（演化）——缺口登记 =========-
-- 归档文本：无向变异提供素材 + 资源环境约束筛选 ⇒ 定向演化。
-- 审计结论：F14 给变异、F13 给竞争，无公理连接二者到"定向演化"。
--   不可导出 → 候选桥接公理 B2（见 GAPS.md）。
--   注意：R5 已证"定向演化"是逐语境命题（Ns≳5 成立 / Ns→0 不成立），
--   B2 若入 Γ 必须带 C_q 资格条件，否则与 R5 冲突。
-- 可证的部分只是素材存在性：
theorem evolution_materials :
    ∀ s : Sys, LifeS s → MatEnergyOK s → CanDivide s →
      (UndirectedHeritableVariation s ∧ Competition s) := by
  intro s hs hm
  exact ⟨F14 s hs (F1 s hs hm).1, (F13 s hs).2⟩

-- ========== 推论 4（耗散结构）——骨架可证，完整版缺口 ==========

/-- 骨架方向 1：生命 ⇒ 捕获能量 ⇒ 废热（F6→F8）。 -/
theorem life_produces_waste_heat :
    ∀ s : Sys, LifeS s → WasteHeat s := by
  intro s hs
  exact F8 s hs (F6 s hs)

/-- 骨架方向 2：有序结构存续必须持续输入自由能（F9 前件）。 -/
theorem structure_needs_energy :
    ∀ s : Sys, LifeS s → NeedsEnergyForStructure s := fun s hs => (F9 s hs).1

-- 完整版（"以环境熵增为代价维持有序"）须 M1 + 生命-环境合系统的孤立性前提；
-- 合系统孤立性不在 Γ_LE 签名内 → 登记为 GAPS.md 缺口 C4-full，不入本文件。

-- ========== 推论 5（生态共演化）——可证 ==========

/-- 双向耦合 ⇒ 环境筛选生命 ∧ 生命改造环境（H(L←E) ∧ H(E←L)）。 -/
theorem coevolution_bidir :
    ∀ l e : Sys, LifeS l → EnvS e → Coupled l e →
      EnvSelectsLife l e ∧ LifeReshapesEnv l e := by
  intro l e hl he hc
  exact ⟨H_LE l e hl he hc, H_EL l e hl he hc⟩

/-- 闭环方向：耦合下生命反向写入环境（E7）。 -/
theorem life_writes_back :
    ∀ l e : Sys, LifeS l → EnvS e → Coupled l e → LifeWritesBackEnv l e :=
  fun l e hl he hc => E7 l e hl he hc

-- ========== 推论 6（吸引子分类）——超出一阶骨架 =========-
-- 归档文本：体系自发抵达四类终态（稳态/振荡/生态锁定/全局灭绝）。
-- 审计结论：终态分类是动力学稳定性陈述，需要状态空间与流结构，
--   FOL 谓词层不可达。处置：交仿真层（x-s-mera 内核的四类终态相图测绘），
--   不在 Γ_LE 增补公理。见 GAPS.md B3。

end XSMERA
