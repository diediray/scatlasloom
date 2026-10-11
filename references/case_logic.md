# 原论文的组织方式与可迁移判断

来源：Nature Communications，DOI [10.1038/s41467-025-59341-w](https://www.nature.com/articles/s41467-025-59341-w)，Fig.2–5及相关结果、Methods “Clustering and cell subtype identification”及“Calculation of gene set activity score”。本次核对使用教程040保存的原始Fig.2、全文XML及作者代码。阅读范围为构图相关结果和方法，并非对整篇研究的全部证据作完整复审；本任务值得精读这些部分，因为它们能直接决定skill如何选列、分层和组织后续问题。

## 先理解这张图在研究中做什么

论文要解释成熟TLS与局部抗肿瘤免疫的关系。Fig.2的起点是：不同TLS病理组别中的T/NK细胞，有没有不同的组成和状态？分群提供比较对象，组别气泡找出关联，状态热图解释这些关联的生物学内容，通路热图补充信号背景。它们共同帮助作者选择后续需要追踪的细胞。

原文结果指出CD8的1/2群具有初始/记忆及stem-like相关特征，在mTLS中富集；4群具有耗竭特征，更关联imTLS/nTLS。仅把这些细胞统称CD8就会丢掉这组差别。Naive、Cytotoxic与耗竭列使比较更有内容，Proliferation则刻画横跨谱系的cycling群。Prog exh列又将CD4第12群与第13群的不同状态区分开，为后续CD4 Texprog/Tfh的研究提供起点。

这五列不是五个互斥细胞类型。细胞毒与耗竭可以共存；作者后续将CD8的3/4群共同用于functional表型分析。因此“耗竭分数高”不能简化成“没有杀伤能力”，“Naive分数高”也不能单凭这一列证明细胞的分化潜能。

证据推进有具体对象：Fig.2d以组织染色检查相关细胞的富集；Fig.2e/f以轨迹和沿轨迹的marker变化研究CD8状态关系；Fig.2g–m加入TCR克隆扩增及共享，支持stem-like与functional群的联系。Fig.3进一步研究CD4 Tex与Treg的轨迹及克隆特征，并关注Texprog/Tfh。组合图本身提供的是亚群关联与表达状态，后续分析各自增加信息。论文也明确承认这些细胞的实际来源及局部生成过程仍有推断限制。

对skill的启示：先问一张图需要区分哪些有意义的解释，才知道哪些列不可少。若合并Prog exh与Term exh，会丢掉本文要追踪的CD4差别；若只显示CD8总体比例，会隐藏内部状态组成的变化。

图的纵向组织是细胞类别及其亚群；横向证据是组别关联、功能状态与信号活动。总UMAP让读者看见全体；灰底小图帮助定位父类；编号与颜色把空间点云和矩阵每一行连接起来。表中同一行必须是同一亚群，不能让各热图分别聚类后失去对应。

原Fig.2a为108,388个细胞、28亚群。以前教程模拟图写80,000，属于样式演示。亚群身份：1–5 CD8 T；6 cycling CD4 T、7 cycling CD8 T、8 cycling NK；9–14 CD4 Tconv；15–21 Treg；22–26 NK/NKT；27 γδT；28 ILC。四个主要高亮面板不是完整、严格的四分分类树。cycling把不同谱系的共同状态单独组织为蓝色块，生物身份与展示分块应各自保存。

机器可读的身份/展示分块对照见[tls_cluster_semantics.tsv](tls_cluster_semantics.tsv)，其中color来自此前发表图取色。高亮集合的独立定义示例见[highlight_membership_examples.tsv](highlight_membership_examples.tsv)：这些是我们解释新接口的设计例子，并非恢复出的作者小图掩码。输入来源与指纹见[source_manifest.json](source_manifest.json)。

群数不用于反推resolution。方法中的总体Leiden为resolution=1；T/NK再聚类用2000 HVG和40个Harmony分量。细胞多、画面复杂并不自动意味着高resolution。

## 各模块确切含义

| 模块 | 内容与作用 |
|---|---|
| TLS组别气泡 | nTLS=无TLS，imTLS=未成熟TLS，mTLS=成熟TLS；病理定义的新鲜肿瘤样本数分别4/4/6。GLM比较组成，颜色log2 OR、面积对应Bonferroni校正后−log10 P。作者代码最终使用eff contrasts；复现前核对实际contrast，不把它随口解释成固定两两比较或合并其余组。 |
| Naive | CCR7、IL7R、TCF7；初始/记忆相关程序。 |
| Cytotoxic | GZMA、GZMB、GZMH、PRF1、NKG7、GNLY；杀伤程序。 |
| Proliferation | STMN1、TUBB、MKI67；增殖程序。 |
| Prog exh | CXCL13、TCF7、PDCD1；祖细胞样耗竭/Tfh相关程序。 |
| Term exh | TOX、HAVCR2、CD101、PDCD1、LAG3、CXCL13、ENTPD1、TIGIT；终末耗竭程序。 |
| 通路 | JAK-STAT、MAPK、NF-κB、PI3K、TNFα、TRAIL、WNT；PROGENy响应基因权重结合MLM推断活性，显示Z-scaled scores。 |

功能热图是五列。图注称平均状态评分，图例写Mean expression；仅凭这两处尚不能确定每一步计算和缩放。精确数值复现需找到对应矩阵/代码。旧教程写成四列、将PI3K误写TGFβ，均已明确纠正。GSE194247参考案例自身可以包含TGFβ，这与纠正原图列名是两件事。

“成熟TLS富集哪些细胞”与“这些细胞呈现何种状态”在一行连接，后续才进入轨迹、克隆共享和组织验证。学习的是这个证据组织方法。该论文Fig.4的B细胞图换成通路及免疫球蛋白表达，Fig.5的髓系图换成通路及MHC表达，说明列的选择本来就随问题变化。

## 原文事实与我们的设计提炼

- 原文事实：nTLS、imTLS、mTLS来自不同样本的病理分组。它们不是对同一患者的纵向时间点；图上的左右排列本身不构成演变轨迹。
- 原文事实：CD8、CD4 Tconv、Treg等身份与cycling等状态在展示上并置，颜色和编号跨面板对应。亚群编号不是拟时序。
- 原图可见：同类亚群常用相关色调，cycling单列蓝色，Tconv红、Treg紫。我们的设计提炼：将类别关系与特殊状态有意识地编码进色系。未发现作者明确陈述“红紫相邻是因为同属CD4”，所以将其记录为设计解读。
- 原文事实：ISG群的JAK-STAT、cycling群的MAPK、Treg部分群的TNFα活动具有特定模式。我们的设计提炼：通路模块要提供与表达状态互补的解释，而不是重复列一组高表达基因。
- 稳定保留：同一亚群的证据对齐、父子层级能辨认、数值定义准确、颜色一致。随问题改变：谱系数、亚群数、对比条件、热图内容、群的展示顺序及是否需要所有模块。

## 用同一论文检验迁移思路

Fig.4研究B细胞成熟及浆细胞生成。Ig重链恒定区表达可帮助比较不同Ig类型，并与成熟状态及后续class-switch相关分析连接；把T细胞耗竭五列照搬过去会错过这个问题。Fig.5的MHC模块则帮助理解髓系亚群的抗原呈递相关特征。二者共享组织方法，但各自具有选列理由。

2026-10-10直接核对三张完整高清原图后补充：Fig.4b的Ig与Fig.5b的MHC均为堆叠小提琴矩阵，颜色图例为Mean expression，虽然图注概称heatmap。Fig.5b通路为九列，较Fig.2/4多TGFβ、VEGF。Fig.5e中cDC2的功能富集进一步衔接到Fig.6的空间与配受体分析。选列、图形及后续结果的对应见[scientific_panel_design.md](scientific_panel_design.md)。

这只是基于原论文的设计推演检查，不等于已经完成独立数据或独立模型的执行验证。

## English summary

The transferable design links a hierarchical cell inventory to condition association, expression programs and signaling activity through common row identities and colors. Fig.2 contains 108,388 cells, 28 subclusters, five state columns and seven pathway columns including PI3K. Cycling CD4, CD8 and NK form a cross-lineage display block; γδT and ILC remain distinct identities. The 80,000-cell tutorial image is synthetic. Fig.3 follows CD4 T-cell states; Fig.4 and Fig.5 change the evidence modules to immunoglobulin and MHC expression. State signatures are not mutually exclusive identities: cytotoxic and exhausted features can coexist. The figure identifies associations worth following with tissue staining, trajectories and clonotype analyses; it does not itself establish lineage transitions or causality. Color kinship is our design interpretation, not a verified statement of author intent. The B-cell and myeloid examples support a design-level transfer exercise, not an independent execution test.
