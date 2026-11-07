# Kim2018 biomarker-stratified response panel (standalone review analysis)

This deterministic tabulation uses the existing frozen patient-level data. Formal manuscripts and original supplemental tables were not modified.

Run: `python analyze_clinical_stratification.py` (defaults resolve from the script location). The saved input snapshot is preferred when present; otherwise the original project data path is used. The manifest separately preserves the original upstream source path and hash. No third-party packages are required.

Verified source totals: 45 patients, 12 responders; MSI and EBV ISH are complete; CPS is available for 41 patients and missing for 4.

| Configuration | Original responders/patients | MSI-positive | EBV ISH-positive | MSS and EBV-negative responders/patients | CPS available | CPS median (Q1-Q3) |
|---|---:|---:|---:|---:|---:|---|
| Dual-low | 3/13 | 1 | 1 | 1/11 | 13/13 | 1 (0-2) |
| Fibroblast-dominant | 1/9 | 2 | 0 | 0/7 | 8/9 | 0.5 (0-1.25) |
| Macrophage-dominant | 8/9 | 3 | 4 | 1/2 | 8/9 | 33.5 (11.5-50) |
| Dual-high | 0/14 | 0 | 0 | 0/14 | 12/14 | 0 (0-0.25) |

The dual-high restricted group still contains 14 patients with zero responders; zero is the response count, not the group size. The macrophage-dominant restricted patients are PB-16-005 (responder) and PB-16-060 (nonresponder).

No new inferential test or response model was fitted. The small restricted group is descriptive and cannot establish biomarker-independent association. No multiplicity family or random seed applies to this deterministic tabulation.

Files:
- `Supplemental_Table6_biomarker_stratified_response_panel.tsv`: publication-friendly metric rows across the four configurations.
- `Kim2018_configuration_biomarker_numeric_summary.tsv`: numeric summary, one row per configuration.
- `Kim2018_response_by_restriction.tsv`: original, MSS-only, EBV-negative-only and joint-restricted counts and response fractions.
- `Kim2018_patients_with_restriction_flags.tsv`: all original columns plus explicit biomarker/restriction flags.
- `Kim2018_patient_input_snapshot.tsv`: unchanged input snapshot.
- `clinical_stratification_manifest.json`: source path, hashes, definitions, runtime and validation.
- `Table6_panel_footnote_EN.txt`: English table footnote.

English table footnote:

Configurations and response labels were taken from the frozen Kim2018 patient-level source data; no signatures or configuration thresholds were re-estimated. MSI status was defined by the pentaplex assay (MSI versus MSS), and EBV status by in situ hybridisation. The restricted analysis included only patients classified as both MSS and EBV ISH-negative. Missing biomarker values were reported separately and were not treated as negative. PD-L1 CPS is summarised among patients with available measurements as the median and first and third quartiles, calculated by linear interpolation (R type 7). Counts and response fractions are descriptive; no additional hypothesis tests were performed. The restricted macrophage-dominant group contained only two patients, including one responder, and does not establish a biomarker-independent association or treatment benefit. CPS, combined positive score; EBV ISH, Epstein-Barr virus in situ hybridisation; MSI, microsatellite instability; MSS, microsatellite stable.
