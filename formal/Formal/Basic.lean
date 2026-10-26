/-!
# Basic：论域与谓词声明（X‑S‑MERA LE 公理体系 Γ_LE 的 L1 编码）

层级声明（v3.3.1）：**L1 对象层**。
公理来源：《X‑S‑MERA 框架：生命–环境同尺度耦合演化完备公理体系（最终归档版）》
——元公理 M1–M2、生命事实公理 F1–F16、环境事实公理 E1–E7、耦合充要条件（第 2 节）、
严格时序（第 6 节）、六条推论（第 7 节，非公理）。

形式化粒度声明（诚实登记）：本编码是 **Horn skeleton**——
公理被编码为谓词间的蕴含关系（正 Horn 子句 + 一条 iff），
物理内容（熵、能量、速率比较、几何结构）被抽象为谓词。
因此本工程检验的是公理集的**逻辑拓扑**（哪些推论可组合导出、哪些缺桥接条款），
不检验物理内容本身。形式化损失在 GAPS.md 逐条登记。
-/

namespace XSMERA

/-- 论域：系统 -/
axiom Sys : Type

/-- 尺度（粗粒化层级标识） -/
axiom Scale : Type

/-- 逻辑时间 τ（投影空间内部推演时间） -/
axiom Tau : Type

/-- 逻辑时间预序（自反、传递） -/
axiom leTau : Tau → Tau → Prop
axiom leTau_refl : ∀ t : Tau, leTau t t
axiom leTau_trans : ∀ {a b c : Tau}, leTau a b → leTau b c → leTau a c

-- ========== 谓词（L1 签名） ==========

/-- 生命子系统 -/
axiom LifeS : Sys → Prop
/-- 环境子系统 -/
axiom EnvS : Sys → Prop
/-- 谱系（F15 对象） -/
axiom LineageS : Sys → Prop

/-- 尺度一致：粗粒化层级一致 ∧ 截断规则一致 ∧ 论域粒度一致（第 2 节） -/
axiom SameScale : Sys → Sys → Prop

/-- 双向稳定耦合（核心定理左端） -/
axiom Coupled : Sys → Sys → Prop

/-- 单向扰动（尺度不一致时的唯一交互形式，第 2 节后果） -/
axiom OneWayPerturb : Sys → Sys → Prop

/-- 孤立系统 / 开放系统 / 受约束（M1/M2 前提） -/
axiom Isolated : Sys → Prop
axiom OpenS : Sys → Prop
axiom Constrained : Sys → Prop

/-- M1 后件：孤立系统总熵不减少 -/
axiom EntropyNonDec : Sys → Prop
/-- M2 后件：局部耗散最优（非全局最优） -/
axiom LocalDissipOpt : Sys → Prop

-- ========== 生命谓词（F1–F16 骨架） ==========

axiom MatEnergyOK : Sys → Prop
axiom CanDivide : Sys → Prop
axiom Heredity : Sys → Prop
axiom BoundedDeviation : Sys → Prop
axiom SelfSimilarGrowth : Sys → Prop
axiom PopExpand : Sys → Prop
axiom LocalCouple : Sys → Prop
axiom SignalExchange : Sys → Prop
axiom CanEmergeHierarchy : Sys → Prop
axiom DamageAccIrrev : Sys → Prop      -- F4 不可逆损伤累积
axiom DamageOverThresh : Sys → Prop    -- 损伤突破阈值
axiom Death : Sys → Prop
axiom CanCaptureEnergy : Sys → Prop    -- F6
axiom EnergyTransportLocal : Sys → Prop  -- F7
axiom WasteHeat : Sys → Prop           -- F8 能量逃逸耗散
axiom NeedsEnergyForStructure : Sys → Prop  -- F9
axiom EnergyCutoff : Sys → Prop
axiom DamageAccel : Sys → Prop
axiom StructureDeterminesEnergy : Sys → Prop  -- F10
axiom EnergyBoundsComplexity : Sys → Prop     -- F11
axiom InteractionLocal : Sys → Prop           -- F12
axiom ResourceFiniteFor : Sys → Prop          -- F13
axiom Competition : Sys → Prop
axiom UndirectedHeritableVariation : Sys → Prop  -- F14
/-- 定向演化（归档推论 3 的结论谓词；新概念申请，见 GAPS.md B2） -/
axiom DirectedEvolution : Sys → Prop
/-- 演化资格（R5 的 C_q：Ns≳5 语境成立 / Ns→0 不成立；"逐语境宏 M"的谓词抽象。
    新概念申请，见 GAPS.md B2） -/
axiom EvoQualified : Sys → Prop
axiom BelowReplacementRate : Sys → Prop     -- F15
axiom LineageDeclineToExtinction : Sys → Prop
axiom PerceivesEnv : Sys → Prop             -- F16
axiom EnvRegulates : Sys → Prop

-- ========== 环境谓词（E1–E7） ==========

axiom ResourceGeneration : Sys → Prop
axiom ResourceSelfDissipation : Sys → Prop
axiom LocalSpatialTopology : Sys → Prop
axiom MemoryHysteresis : Sys → Prop
axiom EndogenousPerturbation : Sys → Prop
axiom DynamicCarryingCapacity : Sys → Prop
axiom LifeWritesBackEnv : Sys → Sys → Prop  -- E7

-- ========== 耦合算子双向作用（第 6 节 H(L←E) / H(E←L)） ==========

axiom EnvSelectsLife : Sys → Sys → Prop
axiom LifeReshapesEnv : Sys → Sys → Prop

end XSMERA
