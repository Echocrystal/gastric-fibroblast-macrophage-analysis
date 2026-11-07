# Figure 1 v2: response annotation QA

This is a new standalone candidate, copied and adapted from the 04 clinical-biomarker figure. It does not overwrite 04, 06 or the formal manuscript.

- Response labels were verified against `best_response`: CR 3 and PR 9 all map to `response_binary=1`; SD 15 and PD 18 all map to 0. No inconsistent or unmapped response labels occur.
- Every MSI and EBV panel retains all 45 patients. MSI-positive has 5 responders among 6; MSS has 7 among 39. EBV ISH-positive has 5 among 5; EBV-negative has 7 among 40. These responder counts are labelled R on the group axes.
- Each CPS panel retains all 41 observed values: 11 responders and 30 non-responders. The 4 missing CPS values comprise 1 responder and 3 non-responders; they are omitted only from CPS panels, explicitly described in the title and legend.
- Response encoding is uniform: orange (`#F1A340`) filled triangles for CR/PR; blue (`#2166AC`) open circles for SD/PD. Responders are drawn last to reduce occlusion. No observations are selected by outcome. Both biomarker-status boxes now use the same pale neutral grey (`#999999` at 0.13 alpha); clinical states are identified by x-axis labels, and colour is reserved for the response points.
- Binary-group jitter is horizontal and display-only. CPS remains at its measured x values on a labelled log1p display axis. All score y coordinates are unchanged. No fitted line or response-derived cutoff was introduced.
- All nine nominal P and BH q values are reproduced and verified against the original 04 statistics snapshot to absolute tolerance 1e-15. The original all-patient six Mann–Whitney comparisons and three Spearman correlations remain the sole tests. No response-stratified or new response association tests were added.
- The title, figure legend and figure footer explicitly distinguish response symbols from the original score–biomarker statistics. All nine comparisons still form one BH family; none reaches q < 0.05.
- Data are traceable point by point in `Kim2018_biomarker_figure_source_data.tsv`, including patient ID, clinical biomarkers, three scores, `best_response`, binary response and response label. Original input, prior statistics, copied-script origin, current script and outputs have SHA256 provenance in the manifest.
- Academic-figure-skill typography/palette/export baselines remain in the copied script. Final vector PDF and SVG plus 300 dpi PNG are saved on an exact 183 × 170 mm canvas. Minimum figure text is 6.5 pt; panel labels, statistical text and external response legend remain separated from observations.
- Rendered PNG inspected: no clipping, no text-data overlap, group counts and response-symbol legend are legible; direct symbols make response encoding interpretable without colour.

Verdict: ready for review as an independent v2 figure. No inferential interpretation has changed.
