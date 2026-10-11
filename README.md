# Single-cell subtype atlas · 单细胞亚群全景图

A Codex skill and R toolkit for publication-oriented single-cell atlas figures: a main UMAP, aligned population highlights and marker evidence, with optional research-driven panels.

用于制作单细胞亚群全景图的工作流与 R 工具。保留清楚的层级、协调的配色、紧凑的矩阵和逐行对应的证据；根据课题选择扩展内容。

## 核心工作流

**核查数据 → 准备候选分群与命名证据 → 用户明确审核 → 设计扩展模块 → 正式成图与检查。**

- **科学审核是前置关卡。** 用户需确认当前分群及命名版本，之后才能计算下游条件关联、功能／通路模块并制作正式组合图。认可美术效果不等于认可科学标签。
- **保留清楚的视觉结构。** 主 UMAP、灰底高亮、marker 和可选模块按群对应；单谱系与父类—亚群层级分别组织。
- **内容可变。** 群数、resolution、算法、基因、功能和模块数量均由数据与问题决定。
- **注释须有证据。** 结合差异 marker、多基因、竞争谱系、样本支持及 QC；证据不足时保留未决状态。
- **结果可回查。** 输出图形、数值定义、输入映射、选择依据和运行代码。

## 安装与调用

将仓库克隆到你使用的 Codex skills 目录中。例如使用默认目录时：

```bash
git clone https://github.com/diediray/single-cell-subtype-atlas.git ~/.codex/skills/single-cell-subtype-atlas
```

如果你设置了其他 Codex 配置目录，将目标路径改为该目录下的 `skills/single-cell-subtype-atlas`。重启或重新加载 Codex 后调用：

```text
使用 $single-cell-subtype-atlas，检查我提供的单细胞数据。
先提交分群和群命名的审核材料，等我确认后，再制作亚群全景图。
```

提供数据位置、研究背景与已有科学确认即可。详细执行规则见 [SKILL.md](SKILL.md)。

## R 工具与输入

R 依赖按使用的入口选择：

```r
install.packages(c("Matrix", "SeuratObject", "ragg", "svglite", "systemfonts"))
```

`grid` 随 R 提供。现有样式使用 Arial；运行环境需要该字体。PDF 导出在 macOS 使用 Quartz，在其他平台需要可用的 Cairo；发布前应检查字体和最终尺寸。

| 工具 | 用途 | 输入说明 |
|---|---|---|
| `prepare_marker_atlas.R` | 从 RNA counts、固定分群/UMAP 和完整差异表准备 marker 证据 | [准备接口](references/marker_preparation.md) |
| `render_modular_atlas.R` | 主 UMAP、逐群小 UMAP—marker 配对及可选扩展 | [组合图接口](references/modular_render_contract.md) |
| `render_atlas.R` | 父类／亚群、独立高亮成员和群级矩阵 | [层级接口](references/render_contract.md) |
| `render_marker_atlas.R` | 按需使用的 marker 与细胞表达热图 | [marker 接口](references/marker_render_contract.md) |
| `score_helpers.R` | 有明确分析设计时使用的表达、通路和配对组成辅助函数 | [数值定义](references/analysis_choices.md) |

紧凑组合图入口目前支持 2–9 群，每群最多六个已选 marker；更复杂的层级使用独立接口或适配实现。它不是自动完成分群和注释的一键分析软件。脚本本身不代替用户的科学审核。

## 合成示例

```bash
Rscript examples/make_synthetic_job.R
Rscript scripts/render_modular_atlas.R examples/generated/figure_job.rds examples/generated synthetic_atlas
```

示例使用生成的坐标和人工特征，仅演示绘图接口，不代表真实 UMAP 分析、细胞身份或生物学结论。

![Synthetic rendering example](examples/synthetic_preview.png)

## 来源与验证范围

证据排列方式参考 [Nature Communications, DOI: 10.1038/s41467-025-59341-w](https://www.nature.com/articles/s41467-025-59341-w)。相关方法与视觉解读保留在 `references/`；仓库不包含论文原图、作者原始脚本、用户的单细胞数据或对话记录。

参见 [验证范围](references/validation.md)。图形工具的运行检查不等于对新数据的注释、统计或投稿合规作保证。
