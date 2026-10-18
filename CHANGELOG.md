# Changelog

本仓库采用语义化版本。理论条款编号（A1–A12、C1/C2、F1/F2、I1）与《X-S-MERA 修订登记》保持同步。

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
