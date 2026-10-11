# 基础marker版接口

调用：`Rscript scripts/render_marker_atlas.R marker_atlas.rds output_dir [prefix]`。

输入由模型从项目准备，为R list：

- title：简短图标题，可省略。
- cells：cell_id、cluster_id、umap_1、umap_2，可加patient_id。
- clusters：有序群ID；colors：以群ID命名的HEX颜色。
- markers：按群及marker次序排列的数据框，group、feature、row_id必需。row_id区分不同群重复出现的同一基因。
- dot_mean_z：marker行×群列，群内mean log RNA按基因跨群标准化。
- dot_pct：相同结构，raw RNA>0的百分比。
- cell_z：同一marker行×全部细胞列，按基因在全细胞log RNA上做z-score；列名与cells的cell_id严格一致。
- metadata：表达层、基因选择、缩放及细胞排序的说明。

主UMAP使用冻结坐标。每个群的小UMAP与该群marker行块居中对齐；点图显示这些marker在所有群的表达；右侧热图同一行是同一基因。全细胞按群、患者、barcode排列，各群热图块等宽，列头写真实细胞数，因此热图块宽不代表丰度。

点图颜色显示±1.5以内的群均值z，点面积表示检出率；细胞热图显示±2以内的细胞表达z，超过范围仅截断显示。平均表达z与细胞表达z的分母不同，各自有图例。

输入保留全部细胞与数值。热图输出为栅格；密集细胞列在有限像素下显示，不应从缩小后的单条像素反推单细胞数量。PNG/PDF/SVG及marker_row_alignment/highlight_alignment记录一同导出。

GSE183904修订版验证了5群、24个marker行、10,196细胞；其中C2为4行，其余各5行。临时一页中不同块采用共同基因行高，保留小图的最小高度。
