# 通用marker准备入口

`Rscript scripts/prepare_marker_atlas.R /path/to/marker_job.R`

模型先检查项目，再生成MARKER_JOB配置；用户不需要填写这些字段。配置保存数据路径、字段映射、marker规则和色板，脚本负责重复计算与输出。

## 输入与适用范围

支持已有RNA counts和固定UMAP的Seurat RDS；也支持项目适配器生成的RDS list，包含稀疏整数`counts`、以barcode为行名的`metadata`及两列`umap`。可以用后一种形式接入其他工具的已完成分析。

完整差异证据表为TSV或TSV.gz，用marker_columns显式对应group、feature、padj及可选auc。feature须与counts行ID一致。若使用Ensembl ID，模型先处理显示名与精确ID映射；不要将某个转录本或多个同名基因的P值任意移给合并后的表达。

脚本重新计算log1p CP10K、群均值、检出率及mean-CP10K log2FC，保留输入差异表的原始列。已有差异表的校正方法写在adjusted_p_definition，不自动都叫BH。没有AUC时按重算效应和检出差排序；有AUC时可按AUC优先排序。

## 配置要点

- object_path、markers_path、output_dir：支持相对配置文件的路径。
- cluster_field、cluster_order、colors：固定分群字段、展示顺序及命名色板。
- assay、reduction：Seurat入口通常为RNA与umap，依据实际对象确认。
- coordinates_path与coordinate_columns：可选冻结坐标表；验证与对象一致后优先使用保存的表值，避免文本序列化的微小差别改变参考图。
- marker_columns：例如`c(group='cluster',feature='gene',padj='p_val_adj')`，或含`auc='auc'`。
- filters：max_adjusted_p、min_pct_in、min_delta_pct（百分比点）、min_log2fc、max_per_cluster，以及可选min_auc/exclude_regex。阈值由具体项目决定。
- patient_field与donor_policy：none/report/filter。没有生物学重复标识时不制造患者支持。有标识时donor_rules规定每侧最少细胞、足够患者数及方向一致性要求。

donor_policy=report只报告支持，不改变候选；filter在有足够可评估患者时依规则筛选，支持不足者明确记为当前数据的描述性候选。某群合格基因少于目标数量时按实际数量保留。若一个都没有，先保存筛选证据，交由模型重新检查设计；不自动放宽所有阈值凑数。

## 输出与衔接

输出input/marker_atlas.rds、cells.tsv、完整筛选证据、选择结果、各群数量、患者支持（适用时）、实际配置和输入指纹。当前已认可纵向骨架优先将bundle放入figure_job的marker字段、按研究目的配置panels，再调用[modular_render_contract.md](modular_render_contract.md)的入口。panels可为空。

以下旧入口仅在确需同基因细胞表达热图时使用：

```bash
Rscript scripts/render_marker_atlas.R /output/input/marker_atlas.rds /output/figures figure_name
```
