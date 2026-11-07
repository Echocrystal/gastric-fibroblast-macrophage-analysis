# Figure and statistics QA

Figure: `Figure_R2_clinical_biomarker_relationships`.

Scientific question: relationships of the three existing continuous scores with MSI, EBV ISH and PD-L1 CPS. This figure does not test clinical predictive value or treatment benefit.

Contract: quantitative grid, three score rows by three clinical-biomarker columns; all nine relationships shown without selection; exact 183 × 170 mm canvas; vector PDF/SVG and 300 dpi PNG.

- Asset routing: the basic-plots directory has a README and no executable simple box/scatter asset. The paired box-scatter production asset and genomic marginal-correlation asset were read and their companion PNGs inspected. Their paired-connection and log-log jointplot semantics do not fit unpaired clinical data, so the plot uses parameter inheritance with explicit classification in the script. No synthetic production data or paired connectors were reused.
- Style: required typography, palette and export baseline blocks are present verbatim. Fonts are Arial/Helvetica/Liberation Sans; minimum figure text is 6.5 pt; top/right spines removed; direct group labels make colour redundant. No legends occlude data.
- Data integrity: 45 unique patients; all score values finite. MSI panels include 39 MSS and 6 MSI; EBV panels include 40 negative and 5 positive; CPS panels include all 41 measured values and explicitly omit only 4 missing values. No other observations are omitted or downsampled. Missing CPS is not set to zero.
- Plot integrity: every patient is shown as a point. Boxes are median/Q1/Q3 with 1.5 IQR whiskers. Binary jitter is deterministic and changes only horizontal display position. CPS uses a log1p display scale with original-unit tick labels; no statistical transformation, fitted line or outcome-based cutoff is used. Shared y limits within each score row support comparisons.
- Inference: patients are the independent units. Six two-sided asymptotic Mann–Whitney tests use tie-corrected variance and continuity correction; three Spearman tests use average ranks for ties and asymptotic two-sided P values. All nine tests form one BH family. Both median differences and rank-biserial effects are retained in the statistics TSV; median differences are positive-minus-negative. The figure reports exact nominal P and adjusted q without significance stars.
- Independent statistical cross-check: base R 4.4.0 `wilcox.test(exact=FALSE, correct=TRUE)` / `cor.test(method="spearman", exact=FALSE)` and `p.adjust(...,"BH")` reproduced all nine P and q values; maximum absolute P difference was 2.220446 × 10⁻¹⁶.
- Export verification: PDF page measured 183.000 × 170.000 mm. PNG is 2161 × 2007 pixels with 299.9994 dpi metadata. PDF text is extractable; SVG retains editable text. Input snapshot, script and output SHA256 values all matched the saved manifest.
- Visual pass: rendered PNG inspected at high resolution; all panels align, no text is clipped, no annotations cover points, small positive groups retain individual observations, and colour contrasts remain distinguishable with direct labels.
- Scientific interpretation: none of the nine comparisons reaches q < 0.05 (minimum q = 0.1739). The nominal sender–CPS relationship is reported with its adjusted q. Weak or non-significant findings are preserved; no signal threshold is used to discard panels or alter data. Small MSI/EBV positive groups limit precision, so non-significance is not evidence of no relationship.

Verdict: ready as a standalone review figure. No formal manuscript or original supplemental output was changed.
