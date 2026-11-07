# 临床常用指标与本研究连续评分的关系（独立新图，未合并正式稿）

本次只比较已有的连续 ecological、sender、receiver score 与 MSI、EBV ISH、PD-L1 CPS 的关系，不按治疗应答选阈值，也不训练预测模型。

MSI/EBV 两列保留全部45例：MSI阳性6例、MSS39例；EBV ISH阳性5例、阴性40例。CPS有记录41例，缺失4例，仅在CPS散点图及相关分析中排除这4例。

9项比较统一进行BH校正。各项q均≥0.1739，没有一项达到q<0.05。

| 连续评分 | MSI阳性−MSS中位差；P；q | EBV阳性−阴性中位差；P；q | CPS Spearman ρ；P；q |
|---|---|---|---|
| Ecological score | -0.148; 0.6283; 0.9424 | +0.119; 0.8425; 0.9601 | -0.170; 0.2893; 0.6510 |
| Sender score | -0.299; 0.8806; 0.9601 | -0.841; 0.0629; 0.2830 | -0.364; 0.0193; 0.1739 |
| Receiver score | -0.118; 0.9601; 0.9601 | +0.784; 0.2405; 0.6510 | +0.086; 0.5923; 0.9424 |

数值上，sender score与CPS呈负相关（ρ = −0.364，nominal P = 0.0193），但9项校正后q = 0.1739；因此只能描述该趋势，不能称为校正后显著关联。Ecological score与三种指标的关联较弱。二分类阳性组只有5–6例，未显著不等于不存在关系。

本图可以展示评分与临床指标的重叠程度，但不能单独回答是否独立预测抗PD-1应答；后者需要另外的、样本量足够的联合模型或外部验证。

复现：使用包含matplotlib/numpy/pandas/scipy的Python执行 `analyze_and_plot_biomarker_relationships.py`。脚本默认优先读取同目录输入快照，所有输出保存于同目录。输入来源、SHA256、具体检验、软件版本与图像尺寸记录于manifest；统计全表为 `Kim2018_score_biomarker_relationships_all9.tsv`。
