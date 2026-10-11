# 已认可纵向骨架：可配置扩展模块

运行：`Rscript scripts/render_modular_atlas.R figure_job.rds output_dir stem`。由模型准备输入；用户只需给数据与研究背景。该入口保留上方主UMAP和下方逐群“小UMAP—本群marker点图”的配对。绘图入口不替模型决定研究方向。

`figure_job.rds`是list：

- `marker`：prepare_marker_atlas.R输出的bundle，含cells、clusters、colors、markers、dot_mean_z和dot_pct。
- `panels`：可为空list，或按从左到右顺序放扩展模块。未提供模块即只绘制基础证据。
- `style`：可选width_mm、hero_mm、band_mm、zlim、marker_point_mm、seed。色标默认±2，始终标注。宽度和密度需看图调整。

## 当前实现的模块

所有模块含type、title、features、width_mm。模块中的cluster顺序必须与marker$clusters完全一致。

**heatmap**：values为feature×cluster的有限值矩阵，当前使用与基础marker相同量纲的逐基因z色标。可用groups=list(标题=列索引)给列加轻量分组标题。不可将OR、原始均值或其它量纲混入此z色标；有不同量纲需求时扩展独立色标实现。

**violin**：values为颜色用的feature×cluster z矩阵；expression为非负的feature×cell表达矩阵，列名和顺序必须与marker$cells$cell_id相同。ylim=c(0,最大表达范围)，必须覆盖全部观测。包含零值，零边界反射KDE，每基因跨群共享带宽，每个轮廓峰宽归一化；图注写清共同表达尺度。表达矩阵缺barcode不能凭长度猜顺序，应由模型核对原来源再补齐。

**replicate_points**：data含cluster、replicate_id、gene、mean；每cluster×replicate×gene最多一项。模型在准备阶段定义独立生物学单位、每单位最低细胞数、条件/区域怎样合并，并保存未筛选和筛选后的表。ylim是全部群/基因共同纵轴范围，unit为量纲，replicate_label为简短图例名（例如Patient、Animal）；没有指定时显示Sample。每点一个重复，中位数横线，不自动作检验、不将缺失重复填零。代码不把患者固定为所有数据的样本单位。

## 已验证范围与选择


当前紧凑实现支持2–9群、每群最多6个已选marker；超出范围或新图形类型时适配几何/实现，再验证，不能删群或删证据来通过接口。图高按群数自动计算，超出247mm会提示检查目标版面。更复杂的父类/亚群高亮可用旧render_atlas.R及其独立成员映射。该实现范围不是科学上的分群或marker数量标准。

PNG/PDF/SVG、行对应表与render_metadata.rds写入指定输出目录；原始数据不改。用统一函数生成各格式，核对标签、同群对齐、颜色/大小含义和缺失重复。通过一次输出检查不等于用户认可新案例的审美或科学内容。
