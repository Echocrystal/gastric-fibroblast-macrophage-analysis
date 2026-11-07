# 07 新图2：独立计数与设计核对

只读取 `02_clinical_stratification/Kim2018_patients_with_restriction_flags.tsv` 的45例固定患者数据；未修改04、06或正式文稿，未新增假设检验或拟合模型。`verify_response_counts.py` 可复现本目录的汇总，来源及输出SHA-256记录在 `verification_QA.json`。

## 数据核对

- 45个唯一患者，12例应答：CR 3例、PR 9例；SD 15例、PD 18例。逐患者核对CR/PR对应response_binary=1、SD/PD对应0，45/45一致。
- MSI状态和EBV-ISH状态完整，无缺失。`configuration_MSI_by_response.tsv` 与 `configuration_EBV_by_response.tsv` 各8行，包括n=0的状态组合；n=0时比例留空，不能读作0%应答。
- `full_and_MSS_EBVnegative_response_Wilson95.tsv` 包含全体及MSS且EBV-ISH阴性限制人群的应答数与两侧Wilson95%CI。区间是各比例的描述性不确定性，不是两组差异或独立预测能力检验。

| Configuration | Full response | MSS / EBV-negative response | Restricted Wilson95%CI |
| --- | --- | --- | --- |
| Dual-low | 3/13 | 1/11 | 0.0162–0.3774 |
| Fibroblast-dominant | 1/9 | 0/7 | 0–0.3543 |
| Macrophage-dominant | 8/9 | 1/2 | 0.0945–0.9055 |
| Dual-high | 0/14 | 0/14 | 0–0.2153 |

联合限制共34例，其中2例应答。CPS缺失4例中1例应答、3例未应答：PB-16-005为Macrophage-dominant/PR，PB-16-041与PB-16-051为Dual-high/PD，PB-16-047为Fibroblast-dominant/SD；4例均为MSS且EBV-ISH阴性。逐患者表与按配置汇总分别见 `CPS_missing_patients.tsv` 和 `CPS_missing_response_summary.tsv`。

## 图形方案建议（5句）

1. 建议用45列患者矩阵，按四个configuration分块，行依次为Response、MSI、EBV-ISH、CPS，固定并标明块内排序规则，使同一患者的应答与生物标志物可直接对齐。
2. Response可保留CR/PR/SD/PD四级并用清楚图例说明CR/PR为responders，CPS使用连续色阶且缺失用独立灰色或斜线，不能把缺失并入0。
3. 旁边配全体与MSS且EBV阴性限制后的response forest，每组显示两种分母、应答数和Wilson95%CI，直接表现限制后样本骤减及区间变宽。
4. CPS个体点图保留每例患者，以形状或轮廓标应答，median/IQR只作描述，并在各配置下单独列缺失数及其应答计数。
5. 这套组合比把MSI/EBV与应答交叉分成四类后做堆积条更可读，因为它保留患者层面的共同出现信息并减少颜色类别，同时用forest明确回答限制人群的应答问题。
