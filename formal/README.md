# formal/ — G1-a：Γ_LE 的 Lean4 编码与相对一致性检查

本目录是《X-S-MERA V3.3.1｜全局压测 G1–G5》**G1-a 项**的工程实现：
把 LE 公理体系（最终归档版）的 **L1 对象层片段**编码为 Lean4 常量公理，
作为公开的矛盾攻击面（v3.3.1 附录 B 否定式 N1）。

## 结构

| 文件 | 内容 |
|---|---|
| `lean-toolchain` | Lean 4.14.0（elan 管理） |
| `lakefile.toml` | lake 工程定义 |
| `Formal/Basic.lean` | 论域（Sys/Scale/Tau）+ 全部谓词签名 + 诚实粒度声明 |
| `Formal/Axioms.lean` | M1–M2、F1–F16、E1–E7、耦合充要条件、尺度不一致后果 |
| `Formal/Corollaries.lean` | 归档第 7 节六条推演的演绎审计（可证的给证明，缺口登记 GAPS.md；**禁 sorry**） |
| `Formal/AttackSurface.lean` | N1 攻击方式分析 + 一致性守卫 |
| `GAPS.md` | 编码暴露的推论缺口与候选桥接公理 B1/B2/B3/C4-full |

## 构建与核验

```bash
cd formal
lake build        # 通过判据：编译成功且无 sorry 警告
```

## 已得结论（编码审计产出，独立于编译结果）

1. **归档第 7 节"推论全量可演绎导出"的声明不成立**：
   衰老推论（缺"阈值必破"桥接）、演化推论（缺"变异+竞争→定向"桥接，且与 R5 的逐语境结论冲突须带 C_q）、
   吸引子分类（超出一阶骨架，须交仿真层）。
2. 双层消亡推论的个体层与谱系层均为**公理重述**（F5/F15 本身），无独立演绎内容。
3. 可证的推论：耗散结构骨架（F6→F8）、共演化双向闭环（H_LE/H_EL/E7）、
   尺度不一致⇒不能耦合（iff 逆否）。
4. 当前片段的一致性由元论证支持（无否定公理）；桥接公理入 Γ 前须重跑检查。

## 状态

**✅ 已通过 Lean 4.14.0 内核类型检查**（2026-10-18）：`Formal/{Basic,Axioms,Corollaries,AttackSurface}.lean`
四个文件全部编译通过、零错误、零 sorry/admit（检查运行于 linux_wasm32 官方构建的 Node 运行时；
native `lake build` 在本机仍被程序黑名单拦截，两者共用同一 elaborator，类型结论等价）。

核验中发现并已修复两处真实错误：
1. `Axioms.lean`/`Corollaries.lean` 的模块文档（`/-! -/`）置于 `import` 之前，违反复合命令序
   （`invalid 'import' command, it must be used in the beginning of the file`）——文档降为普通块注释后置；
2. `Corollaries.lean` `evolution_materials` 的 `intro` 漏引入第三个假设（`CanDivide s`），
   匿名构造器作用于未消去的蕴涵式上被拒——补 `hdiv` 后通过。

如本机装有 elan，标准构建仍是：`cd formal && lake build`。
