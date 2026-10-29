# B3-ext · LE 内核四类终态机制扩展

> **状态**：✅ 已实现（v1.6.47）。B3（v1.6.45）Part A 发现当前 LE 内核
> （v3.4-le.1）单调稳定、四类终态不可达 → 登记 B3-ext；本轮为内核引入**可选
> 扩展机制**，使四类终态在 **LE 内核本身**可达（非扩展原型）。

## 动机（承接 B3 发现）

B3 Part A：`le_coupling_kernel` v3.4-le.1 耦合网格 36 点**全部收敛唯一非零稳态**——
内核是单调稳定系统（逻辑斯谛增长 + Michaelis-Menten 摄食 + 线性资源供给），
**振荡/灭绝/锁定不可达**。Part B 仅以实验内扩展原型证明四类可达（非内核本身）。
⟹ 登记 **B3-ext 待办**：须在内核引入机制。

## 本轮：可选机制（默认 None 保持向后兼容）

`src/le_coupling_kernel.py` 新增 `mechanisms` 参数与 `EXT_MECH_KEYS` 校验。

| 终态 | 机制键 | 动力学原理 |
|---|---|---|
| **全局灭绝** | `death_threshold`/`death_rate` | life 低于死亡阈值 → 不可逆加速死亡（吸收态） |
| **生态锁定** | `allee_threshold`/`allee_strength` | Allee 双稳：bio 低于阈值趋 0、高于趋 k → 低态/高态锁定 |
| **周期振荡** | `oscillation_gain` + `env_r`/`consumer_mortality`/`half_override`/`uptake_override`/`conversion` | **纯消费者 R-M 结构**：bio 完全依赖资源摄食（`bio'=e·a·bio·R/(h+R)−m·bio`）+ 资源 logistic 再生 → Hopf 极限环 |
| **稳态** | （默认，无机制） | 原单调稳定内核，行为不变 |

## 实验结果（`b3_ext_mechanisms.py` → `results.json`，seed=20261029）

**四类代表用例全部可达**：

| 代表 | final_bio | 终态 |
|---|---|---|
| steady（无机制） | 1.309 | 稳态 |
| extinction（death_threshold=0.90） | 0.000 | 全局灭绝 |
| lock（allee_threshold=1.00） | 0.164 | 生态锁定 |
| oscillation（conversion=0.8） | 8.167 | 周期振荡 |

**机制分岔扫描**（操作化相图）：

- **灭绝阈值分岔**：death_threshold <0.9 → 稳态；≥0.9 → 灭绝。
- **锁定双稳分岔**：allee_threshold ≤0.9 → 稳态；≥1.0 → 锁定（低态）。
- **振荡 Hopf 分岔**：conversion 0.5 → 稳态；≥0.8 → 持续振荡（极限环）。

## 验证与诚实边界

- `python -m pytest tests/test_b3_ext_mechanisms.py` **9 passed**
  （默认向后兼容稳态 + 四类可达 + 三处分岔 + 非法机制键拒绝）。
- 全套回归见 CHANGELOG v1.6.47（须确认内核改动无回归）。
- **诚实声明**：机制为**操作化原型扩展**，非对现实的生态学断言；四类标签是对
  数值轨迹的操作化分类；**默认内核（无机制）行为保持 v3.4-le.1 不变**（向后兼容）。
  振荡（Hopf）是纯消费者 R-M 结构的数学性质，不代表认知系统必然振荡。

## 运行

```bash
python b3_ext_mechanisms.py            # 测绘 → results.json
python -m pytest ../.. -k b3_ext       # 测试（9 项）
```
