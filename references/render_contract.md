# 灵活渲染器输入约定

运行入口：`Rscript scripts/render_atlas.R /absolute/path/config.R`。输入整理和配置由模型承担。

`cells.tsv`：cell_id、cluster_id、umap_1、umap_2。`clusters.tsv`：cluster_id、label、color，可加group_id/group_label作为连续展示分块，另存lineage/state/full_label保留身份。行序取clusters.tsv原顺序。可选label_positions.tsv的cluster_id/x/y仅指定文字位置。

小图高亮的成员映射独立于group_id。设置`mini_mode="custom"`及`highlights_file="highlights.tsv"`，表中每行包含highlight_id、label、cluster_id，可加anchor_group。一个cluster可在不同意义的小图中出现，例如“全部CD8”与“所有cycling”，每个小图内部一个细胞只画一次；同一highlight内重复成员或未知cluster会报错。

条带模式按highlight首次出现的顺序排列。侧栏模式还需要anchor_group指定小图对齐到哪个矩阵块；它仅影响位置，不决定成员。不同小图使用不同的有效anchor，支持输入顺序与矩阵块顺序不同。可以只给部分矩阵块配置侧栏小图。

```text
highlight_id  label               cluster_id  anchor_group
cd8_all       CD8 incl. cycling   1           CD8
cd8_all       CD8 incl. cycling   2           CD8
...           ...                ...         ...
cd8_all       CD8 incl. cycling   7           CD8
cycling       Cycling T / NK     6           Cycling
cycling       Cycling T / NK     7           Cycling
cycling       Cycling T / NK     8           Cycling
```

上例只节选几行说明含义，实际文件是完整的制表符分隔表。各highlight内label和anchor保持一致。父类和状态的解释分别写进标题。输出highlight_membership.tsv和highlight_summary.tsv记录真正用于点云的成员及唯一细胞数。

模块可以任意组合。热图长表为cluster_id/feature/value；气泡长表为cluster_id/feature/color_value/size_value。每个模块包含配置中所有cluster×feature且各一次；缺失值需要先查明来源。

```r
ATLAS <- list(
  title="Dataset · cell population",
  subtitle="Research question or provenance",
  prefix="atlas", input_dir="input", output_dir="output",
  layouts=c("portrait","wide"),
  mini_mode="cluster", # cluster / group / custom / none
  mini_placement="strip", # strip / sidebar (group or custom)
  modules=list(
    list(type="bubble", file="markers.tsv", title="Identity markers",
         features=c("GENE1","GENE2"), labels=c("GENE1","GENE2"),
         limits=c(-1.5,1.5), ticks=c("≤−1.5","0","≥1.5"),
         legend="Mean expression (z)", ramp="or",
         size_max=100, size_breaks=c(25,50,100), size_legend="Cells (%)",
         semantics="Mean log1p CP10K z-scored across clusters; percent raw UMI-positive cells"),
    list(type="heatmap", file="programs.tsv", title="Expression programs",
         features=c("program_a","program_b"), labels=c("Program A","Program B"),
         limits=c(0,1), ticks=c("0","0.5","1"), legend="Relative mean", ramp="function",
         semantics="Mean panel expression, min-max across clusters per panel")
  )
)
```

示例仅说明接口，基因和程序由实际任务决定。ramp可用or/function/pathway，分别对应既有的三条色带。额外统计气泡作为另一个模块，独立定义颜色、大小与单位。父类高亮使用所属子群的原有颜色。

`main_mm`可调整主UMAP显示尺寸（默认110 mm），`particle_factor`可覆盖默认密度倍率。此次28群结构示例采用145 mm及0.65，保留输入坐标，改善原先过淡的点云。尺寸参数是可读性选择，不能取代看图。

渲染器依赖R的grid、ragg、svglite、systemfonts以及配置字体。输出PNG/PDF/SVG和行坐标、标签位置、成员映射、运行记录。28群/80,000模拟细胞的纵横侧栏版与纵向条带版已实际运行并检查；尚未据此宣布新真实数据验证完成。具体范围见[validation.md](validation.md)。

## 可选密度与网页预览参数

`row_height_mm`（默认5.8，至少4）、`group_gap_mm`（默认7，至少0）、`sidebar_min_mm`（默认34，至少25）与 `legend_space_mm`（有气泡默认43，至少32；纯热图默认20，至少20）控制展示密度。压缩后需检查小图标题与相邻分块是否重叠；不得靠删群适应版面。

`preview_dpi` 可设为72–450，在标准PNG/PDF/SVG之外输出同一绘制函数生成的 `_preview.png`，用于网页展示。省略上述参数时保留原有布局与导出行为。可复现例子见 [层级演示脚本](../examples/make_hierarchy_demo.R)。
