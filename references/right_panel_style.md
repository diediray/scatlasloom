# 原案例右侧：点图、热图和小提琴矩阵的美术组织

依据：DOI 10.1038/s41467-025-59341-w，已逐一查看完整Fig.2/4/5及b面板；作者公开代码 `figure_code/Fig2_T_NK_dotplot_odds ratio_for_cell_cluster_type.R`。以下分清原图可见特征、栅格近似测量及我们复用时的选择，不能把外观估计写成作者的原始参数。

## 共同的排版基础

三个模块共享亚群顺序和行带。亚群名称右对齐，右侧紧跟固定大小的身份色圆点，再进入条件气泡列。颜色圆点属于身份标记，大小不随统计量改变。状态/通路矩阵与对应气泡共用行中心；小提琴在同一行带内使用靠下的基线，不能机械地把每个轮廓的几何中心移到行中心。

背景白，无外框、网格线或交替行底色。主要依靠直线对齐、清楚的列间距和模块间留白组织结构。组块内部紧凑、组块之间有呼吸空间：Fig.2按类别切开热图；Fig.4/5则连续排列行，不应强制复制Fig.2的所有断块。

列标题放在上方，约45°倾斜、向右上方展开；各模块分别占自己的标题空间，不叠加笨重的标题框。行名与列名以常规字重的无衬线字体为主，面板字母用粗体。未从原始矢量文件确认字体；复用可沿用现有Arial，但不能声称原作者精确字体/字号已核实。

模块下方分别放短横向色条，名称在上、少量刻度在下。点大小图例用中性灰圆，避免与效应方向的颜色混淆。每个模块宽度由列数和信息形式决定，不等宽硬分画布。

## 左侧条件点图

- 三列条件间距均匀，表头下方有薄的分类色条，起到定位作用；正文气泡颜色表达统计效应，不沿用条件类别色。
- 气泡为实心圆，无醒目的黑描边。大点与邻行、邻列都保持空隙，小点仍形成节奏；没有箱框和坐标轴干扰。
- 最大圆直径在Fig.4/5约为行距的0.60–0.65。这个比例可作为新图试排起点，最终取决于列距、大小变量分布和输出尺寸。
- 原公开代码用 `rev(brewer.pal(9, "RdBu"))`、`geom_point()`、`scale_size_continuous()`，并去掉轴线、刻度、网格、背景和边框。该代码自身的方向与图例摆放不完全等同最终拼版，不能直接认定运行该脚本即可得到原文成图。
- 保留大小图例和色条的独立含义。作者代码设置了最小点且对数值封顶；复用时写清大小映射、最小点及封顶，不把“用了size变量”直接写成从零开始的严格面积正比。

## 中间热图

- 每个单元格是无描边的纯色矩形，相邻格直接接合。Fig.4/5的格子接近方形，宽约为行距的1.03–1.06倍；Fig.2也保持相近比例。避免少量行配巨大宽画布形成长条格。
- 大多数单元格使用浅蓝、米白、浅杏等低饱和中间色，少数强信号落在深蓝或酒红，形成视觉重点。颜色由实际数值决定，不能为了柔和而选择性压低真实强信号。
- 原图同时使用了不同渐变：条件效应较柔和，通路色阶深端更深。状态评分/表达量与通路的量纲和色标位置分别保留，不能只因同为蓝红就共用一个归一化。
- 不逐格写数字、不重复行名，不添加自动聚类树来破坏横向对应。少量必要的数值由图例交代。

## 右侧小提琴矩阵

- Fig.4的Ig、Fig.5的MHC是逐群、逐基因的小提琴阵列。每一列是基因，每一行是亚群；一个小提琴内部使用同一填充色，色条标注Mean expression。
- 轮廓纤细、以填充为主，没有粗黑描边、内部箱线、中位数点或细胞散点叠加。低表达的形状可呈窄针或贴近底部的扁形，较宽的形状与邻列留出白缝。
- 纵向表达分布与横向密度轮廓被压缩到各自的行带，基线与可用高度保持规则。低信号处的水平扁形可能来自密度轮廓，不应擅自当作附加轴线逐格描线。
- 当前未核实原作者的核密度带宽、截断、最大宽度归一化和纵向尺度规则。实现时明确这些设置，保持可比较的尺度；不能按每格随意拉伸到同样大小，以“填满格子”替代表达分布。
- Fig.5在MHC-I/MHC-II上方加细横线，再区分classical/non-classical；它是轻量的列层级标注，不是加粗边框。只有数据确有列层级时复用。

## 颜色测量发现

对Fig.4/5图例内部进行栅格取样得到：条件效应与通路色条的近白段约在全长48–52%；表达量色条的近白段约在11–15%，中段已经是杏橙色。这里描述的是已发表图的外观，不说明作者采用了哪一种变换。

Fig.4色条内部0/25/50/75/100%位置的近似颜色：

| 色条 | 0% | 25% | 50% | 75% | 100% |
|---|---|---|---|---|---|
| 条件效应 | #4378B3 | #9FCEE0 | #F6F6F5 | #E8A485 | #B92D35 |
| 通路 | #224273 | #82BBD4 | #F6F5F3 | #DD8068 | #6D0718 |
| 表达量 | #3969A7 | #F3E4D6 | #E5987A | #C43C3F | #720E21 |

这些为出版栅格色条内部的取样值，受到边界、抗锯齿等影响，不是恢复的原始HEX定义或全部端点。新图的量程及中性值必须由数值含义决定，不盲抄白色位置。身份色、条件色和连续色阶各司其职。

## 下一轮试排的验收重点

先用已确定的真实输入试排三个模块的局部，保持已有科学内容。看亚群行是否对应、气泡是否拥挤、热图格子是否紧凑、小提琴是否清楚且有间隙、表头和图例是否形成统一层级。再以最终输出大小查看整图；既看放大后的轮廓，也看缩小后的阅读节奏。完成样式试排不等于已决定新增科学分析。

## English summary

The shared unit is a cluster row band: dots and tiles use row centers, while violins sit above a consistent lower baseline within the same band. The reference uses solid borderless dots, near-square adjoining heatmap tiles, and compact filled violin matrices without heavy outlines, boxplots or point overlays. White space and alignment carry structure; column headers are angled and each module has a compact horizontal legend. Effect, pathway and expression ramps differ, with the expression ramp's near-white region close to its lower end. Pixel measurements describe published appearance, not recovered source parameters; violin smoothing and normalization remain unverified. Preserve scientific scale definitions while borrowing proportions and visual hierarchy.
