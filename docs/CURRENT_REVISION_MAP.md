# Current revision entry points (7 November 2026)

Run commands from the repository root. Keep the relative directory layout.
These commands operate on the supplied snapshots and write inside their module
directories. Use a working copy when comparing with frozen delivered results.
Required Python packages for the listed numerical/plotting modules are numpy,
pandas, scipy, matplotlib, Pillow and h5py. These are package names, not a
validated environment lockfile. Original R/legacy dependencies are recorded
in the original scripts and documentation.

| Result | Entry point under revision_20261103/ | Input and scope |
|---|---|---|
| Composite / receiver / ITGB2 comparison | 01_spatial_endpoint/run_spatial_endpoint_comparison.py | Frozen section/donor coefficients; original cohort analyses |
| Clinical configurations and restriction | 02_clinical_stratification/ | Patient snapshot, descriptive summaries and response counts |
| Programme score / GSEA discrepancy sensitivity | 03_phagocytosis/run_phagocytosis_analysis.py --snapshots-only | Cached public-data-derived union-gene counts and saved metadata; four score tests and original/matched 14-set GSEA |
| Clinical biomarker tests | 04_clinical_biomarker_relationships/ | Nine predefined score–biomarker tests |
| GSEA running curves | 11_Figure4D_four_pathway_GSEA_curves_candidate_20261104/plot_Figure4D_four_pathway_GSEA_curves.py | Frozen 8,187-gene ranking, original GSEA statistics, fixed gene sets; does not rerun inference |
| Overall spatial pooling | 12_SPON2_ITGB2_independent_cohort_pooling_REVIEW_20261104/run_independent_cohort_pooling.py | Reads the same 01_spatial_endpoint inputs for integrity checks; 19 + 9 participants; ST_R excluded from Overall |
| Current S10D | 13_insertion_ready_figures_20261104/spatial/draw_spatial_insertion_candidate.py | Original cohort and new participant-equal Overall results from module 12 |
| Current S14 | 13_insertion_ready_figures_20261104/clinical/plot_S14_insertion_ready.py | Colocated source_data/figure2; response, MSI/EBV status, CPS, full vs restricted population |
| Current S15 | 13_insertion_ready_figures_20261104/clinical/plot_S15_insertion_ready.py | Colocated source_data/figure1; score–MSI/EBV/CPS associations and response symbols |
| Current S5 A–D | 14_FigureS5_GSEA_four_panel_REVIEW_20261104/plot_FigureS5_GSEA_four_panel.py | Colocated curves, cross-checked against module 11/source_data; four-panel display only |

Example commands:

```sh
python revision_20261103/03_phagocytosis/run_phagocytosis_analysis.py --snapshots-only
python revision_20261103/11_Figure4D_four_pathway_GSEA_curves_candidate_20261104/plot_Figure4D_four_pathway_GSEA_curves.py
python revision_20261103/12_SPON2_ITGB2_independent_cohort_pooling_REVIEW_20261104/run_independent_cohort_pooling.py
python revision_20261103/13_insertion_ready_figures_20261104/spatial/draw_spatial_insertion_candidate.py
python revision_20261103/13_insertion_ready_figures_20261104/clinical/plot_S14_insertion_ready.py
python revision_20261103/13_insertion_ready_figures_20261104/clinical/plot_S15_insertion_ready.py
python revision_20261103/14_FigureS5_GSEA_four_panel_REVIEW_20261104/plot_FigureS5_GSEA_four_panel.py
```

For the signature-overlap recount, use its own module README; it uses Python's
standard library. It preserves all seven final comparisons and shared genes.

## Archived sources and limits

The 103 original scripts in code/ are preserved byte-for-byte from the staged
release, with their historical path settings. Their integrity and syntax do
not prove a fresh raw-data rerun. See RELEASE_READINESS.md for upstream gaps.
The 20261103 README and release documents describe the earlier archive; the
current README and this map supersede the earlier statement that no new
phagocytosis/GSEA revisions are included.

`code/figure_assembly/compose_FigureS10_ABCD_20261104_original.py` preserves the
exact cosmetic assembly source. Its paths (project root, font and pdftoppm)
need configuration before execution elsewhere. The corresponding source PDFs
are included under inputs/; the original A–C filename was FigureS11, which is
historical numbering. Statistical source code is separate from this PDF
composition step. No new scientific calculation occurs during composition.
