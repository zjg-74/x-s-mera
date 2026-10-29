# A13 实例申报表 —— 约定声明字段模板（OPEN-3 落地）

> **状态**：✅ 已落地（v1.6.48）。承接 OPEN-3 连续化规格评审定稿（v1.6.31，
> `docs/xs-framework/32-v3.4-open3-continuous-observability-spec.md` §6）的
> 实现路线第 3 条：**把 θ_min/Δ 纳入实例申报表约定声明字段模板（A13 申报规范）**。

## 0. 一句话

OPEN-3 把 B1 维"开/闭二分理想化"连续化为可观测度 `θ_k ∈ [0,1]`，可读判据连续化为
`∃k: θ_k ≥ θ_min`。θ_min 与 Δ 是**观测协议声明的约定性输入**（与 OPEN-2 A/B、z、η
同类），框架不裁定取值、只要求显式申报。本文件为 A13 实例申报表补齐这些字段。

## 1. 约定声明字段模板（与 OPEN-2/z/η 并列登记）

| 字段 | 键名 | 取值/约束 | 语义 | 同类先例 |
|---|---|---|---|---|
| 约定 A/B | `open2_declaration` | `"A"` / `"B"` | C 级自足读的约定：A=工具算界内自足为定义上的真；B=相对于观测协议的条件性结论 | OPEN-2 |
| z 顶点带阈值 | `z_vertex` | 正数（约定 1.96） | R34 判读禁区 `\|β̂\| ≤ z·SE(β̂)` | OPEN-1/R34 |
| η 尺度阈值 | `eta` | `[0,1]` | R-043 落带 + R44 可靠子区间取值原则 | R-042/043/044 |
| **θ_min** | `open3_theta_min` | `[0,1]` | C 级辨识层可读阈值：`θ_k ≥ θ_min` → 可读 | OPEN-3 |
| **Δ 不确定带** | `open3_delta` | `[0, θ_min)` | 边缘不判带宽度；`Δ=0` 退化回二值 | OPEN-3 |

**约束**：`open3_delta ∈ [0, open3_theta_min)`；`open3_delta = 0` 等价于原二值开闭。

## 2. OPEN-3 连续化 C 级辨识层三值判读（R34 先例降级）

| θ_k 区间 | 判读 | C 级产出 |
|---|---|---|
| `θ_k ≥ θ_min` | **可读** | 进 C 级辨识层，可补自由度 |
| `θ_k ∈ [θ_min − Δ, θ_min)` | **边缘可观测，不判** | 登记"待更强观测"，不产出 C 级闭合结论 |
| `θ_k < θ_min − Δ` | **闭** | — |

**诚实边界**：Δ 取值为约定输入，框架不裁定；边缘带"不判"不产出任何闭合结论。
（若弱可观测导致 E/U 序列噪声增大，属 OPEN-4 噪声方差线，不属本模板。）

## 3. 关键定理（OPEN-3 连续化不改变分级）

- 连续化只**加一个约定阈值参数 θ_min**（+ Δ 边缘带），不新增 L1 公理、不新增判定逻辑层。
- **A/B 级分级不受影响**【定理】：连续化是 C 级辨识层的参数化，B 级判读（z 容差带）不变。

## 4. A13 示范填充（见 `a13_contrast.json` 顶层 `declaration`）

A13 两个场景（A_drill / B_hollow）共享一组约定声明，示范字段取值：

```json
"declaration": {
  "open2_declaration": "A",
  "z_vertex": 1.96,
  "eta": 0.05,
  "open3_theta_min": 0.8,
  "open3_delta": 0.1,
  "c_level_rule": "theta>=theta_min readable; theta in [theta_min-delta,theta_min) edge-nopad; theta<theta_min-delta closed"
}
```

> 该块为**示范性约定声明**（演示如何申报，不构成对现实的测量断言）；θ_min/Δ 取值
> 由观测协议声明，框架不裁定。

## 5. 复现

```bash
# a13_contrast.json 已含 declaration 块；校验：
python -c "import json; d=json.load(open('experiments/a13_first_case/a13_contrast.json')); assert 'declaration' in d; print('declaration present:', d['declaration'])"
```
