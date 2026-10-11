---
name: single-cell-subtype-atlas
description: Create scientifically informative single-cell subtype figures with a main UMAP, aligned cluster highlights, systematically selected markers and research-question-driven evidence panels. Use for 单细胞亚群全景图, UMAP配各个亚群marker, or the tutorial 02/040 style. Preserve both biological information gain and polished visual organization.
---

# 单细胞亚群全景图

用清楚的亚群证据回答具体的生物学问题，并保持参考案例的视觉组织。基础marker说明各群由哪些基因区分；右侧进一步呈现有信息增量的比较或功能特征。用户提供数据路径和已有背景即可，模型负责阅读背景、选择证据、整理输入和成图。

## 必须先由用户审核分群与命名

**分群结果与群命名必须交给用户仔细审核，获得对当前版本的明确确认后，才能继续下游分析和正式组合图。** 这是用户明确要求的科学审核关卡；美术、配色或构图得到认可，不等于分群或命名通过。要求“继续分群”“按案例做图”也不能代替对随后产生的候选结果的确认。

审核前可以完成数据核查、候选分群及参数对照、用于注释的差异marker/共表达/样本支持，并准备审阅图。向用户提交可直接核查的结果：当前分群版本与细胞归属、候选名称和父类、支持与冲突marker、样本/QC情况、拟合并或排除项及未决问题。随后停在科学审核阶段，不能先计算右侧条件富集、功能/通路等下游模块或完成正式组图，再补问用户。

在项目中分别记录分群、命名和视觉认可状态，并关联具体的细胞归属表、名称表版本以及用户确认。已有明确科学批准且输入未变时，沿用该批准，不重复询问。重新聚类、合并/拆分、改变纳入细胞、父类归属或群名后，对变更部分重新提交用户审核；仅调整颜色、字号、位置和留白无需重新审核分群。科学审核未完成时保持待审核，不因等待时间或视觉认可而自动推进。

## 固定基础骨架，按课题确定右侧内容

用户已认可的基础骨架是上方主UMAP、下方B图中逐群的小UMAP紧邻本群marker点图。保留这一配对、对齐与阅读顺序；固定的是证据组件，marker名单与合格数量仍由新数据决定。

右侧是按研究问题设计的扩展区。基因、功能分类、比较对象、模块数量和图形形式均可变化；不将某次案例的热图＋小提琴固定为必备组合。具体流程见[scientific_panel_design.md](references/scientific_panel_design.md)：背景清楚时先向用户简要汇报拟展示的内容及科学目的，再推进；背景不足、不能判断后续要研究什么时，先询问用户，得到方向后再决定扩展内容。等待期间可继续基础marker工作。

## 先完成基础证据

读取项目说明、固定分群/UMAP、物种、完整cluster-vs-rest差异表与raw RNA对象。保留已认可的细胞、编号和坐标。若没有完整marker表，先用适合现有数据的差异分析生成，再组织展示。

优先读用户项目的单细胞教程和既有marker选择规范。每群独立筛选候选，综合正向效应、群内检出率、群间检出率差、AUC（如有）、多重校正结果和样本支持。以3–5个/群作为紧凑展示的起点；证据不足时保留实际合格数量并解释，不按固定总数随意挑几个基因。完整规则见 [references/marker_selection.md](references/marker_selection.md)。

保存候选及最终marker表，记录阈值、排名、排除理由、数值定义。技术基因或应激信号怎样处理须服务于本次目的：展示时可排除主要反映线粒体/核糖体的条目；确实区分当前群的应激、免疫或其他谱系信号应如实呈现。差异marker描述当前群，进一步命名需要组织背景和多基因证据。

可使用 `scripts/prepare_marker_atlas.R` 完成可重现的准备。模型按项目生成配置，映射字段并明确筛选/患者支持规则；支持Seurat或稀疏矩阵适配对象。输入范围和配置说明见 [references/marker_preparation.md](references/marker_preparation.md)。

项目关注的FAP等基因可以另作有明确目的的补充。每个群的代表marker仍由完整差异证据决定。

## 小UMAP后面紧跟该群的证据

默认采用以下读图顺序：主UMAP概览 → 左侧逐群高亮 → 同一水平块内的marker及其证据。不要把小图全部排在远离相应矩阵行的底部。

对于单层分群，推荐每个群占一个水平块，内部是该群的marker行；点图列为全部群，检查marker相对特异性。单细胞表达热图适用于需要辨认群内分布或共表达的情况；若只重复点图的主要信息，则不占用主图右侧。旧marker渲染器会导出该热图，不能因此把它纳入新图。左侧高亮图与该群marker块居中对齐。这里每行代表基因，明确标注，不能与“每行一个亚群”的状态矩阵混淆。

对于有父类/亚群层级的数据，可按父类高亮，对应多个亚群行。谱系身份、跨谱系状态及展示分块分别保存。高亮成员可单独定义，比如cycling矩阵块与全部CD8高亮不必有相同成员。

当前已认可的纵向骨架使用 `scripts/render_modular_atlas.R`，输入见 [references/modular_render_contract.md](references/modular_render_contract.md)。基础marker与可选扩展模块分开传入；panels为空可只画基础证据。当前支持逐基因z热图、细胞表达小提琴、生物学重复点图；新图形按研究问题扩展，不把当前实现当成分析菜单。需要细胞表达热图时可参考旧 `render_marker_atlas.R` 及 [marker_render_contract.md](references/marker_render_contract.md)；复杂父类/亚群高亮使用 `render_atlas.R` 和 [render_contract.md](references/render_contract.md)。

## 从具体问题决定扩展

本节在当前分群与命名已经通过上述用户审核后执行。右侧方向清楚或已经获准，不能跳过科学审核关卡。

基础展示先回答“有哪些群、它们由哪些基因区分”。条件富集、功能基因集和通路是针对进一步问题的可选分析，不按固定四模块填满画布。原论文的TLS、耗竭及PROGENy选择有其研究问题；新数据不自动继承这些列。

需要扩展时，先读 [references/scientific_panel_design.md](references/scientific_panel_design.md)，用简短figure_design.md写清结果问题、各模块增加的具体信息及与下一步分析的联系，再选择方法和证据。先写待检验问题，计算后据真实结果改写结论。科学价值可以是分清状态、比较条件关联、解释功能差异或筛出值得追踪的对象；不要求每块都建立机制。组成差异、marker表达、状态评分和通路活性分别定义单位。只因元数据里有条件，或因为模型文件已经可用，不足以决定新增分析。

原图科学逻辑见 [references/case_logic.md](references/case_logic.md)，跨T/NK、B及髓系的设计例子见 [references/design_examples.md](references/design_examples.md)。需要统计/评分时再读 [references/analysis_choices.md](references/analysis_choices.md)；`score_helpers.R`为可选数值工具。

## 精修与交付

类别色按ID贯穿全部面板；同类或相关类别用协调色系，群内保持可辨。保留细点、灰底高亮、等比例坐标、白色标签描边、统一行距和清楚的图例。样式参数在 `assets/style/`，细节见 [references/visual_grammar.md](references/visual_grammar.md)。

直接查看PNG及PDF，检查小UMAP与证据块对应、基因顺序、点大小与颜色、热图缩放、可读性和裁切。同时将主图与原案例并排查看：行对应、点云质感、分块留白、矩阵紧凑程度和文字层级均须清楚。科学内容和视觉质量分别检查，不能互相替代；优先减少冗余模块，避免靠缩小文字容纳更多内容。

交付主图、PNG/PDF/SVG、marker选择表、准备/绘图脚本与简短说明。报告基因选择和真实观察，不把程序或通路强行解释成新身份。任务反馈不符合目的时修正设计并更新默认入口，保留旧数值作历史记录。

本 skill 借鉴公开论文中的证据组织与本地开发经验。实现范围和验证边界见 [references/validation.md](references/validation.md)。
