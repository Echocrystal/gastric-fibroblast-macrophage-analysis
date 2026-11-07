> Updated 7 November 2026: CURRENT_REVISION_MAP.md identifies the current final plotting entry points.

# Revision modules in this source archive

Use Python 3.12 with NumPy, pandas, SciPy, matplotlib and threadpoolctl for the
scientific plotting modules; the clinical-count table script uses only the
standard library. Installed versions at analysis time are recorded in the
per-module manifests. These are provenance records, not a fully validated
global dependency lockfile.

From the repository root, example commands are:

```sh
python revision_20261103/01_spatial_endpoint/run_spatial_endpoint_comparison.py
python revision_20261103/02_clinical_stratification/analyze_clinical_stratification.py
python revision_20261103/04_clinical_biomarker_relationships/analyze_and_plot_biomarker_relationships.py
python revision_20261103/07_response_annotated_figures_20261103/figure1/plot_biomarker_relationships_response_v2.py
python revision_20261103/07_response_annotated_figures_20261103/figure2/draw_response_integrated_configuration_v2.py
```

Outputs are written within the module folder unless an explicit output option
is supplied. Work in a copy when preserving the archived numeric outputs.
The spatial module needs the retained frozen `inputs/` and reference
`scripts/` files; it does not require source-project paths when those files
already exist. The figure2 script depends on the preserved sibling 02 and 04
directories, so keep the revision hierarchy intact.

The response-annotated clinical figures are candidate revision displays;
response labels add descriptive patient information. The nine score–biomarker
P/q values retain their original all-patient test family. Full and restricted
clinical populations overlap; no new response-prediction model was fitted.
The current phagocytosis sensitivity, GSEA and final plotting modules are
now included. See CURRENT_REVISION_MAP.md for their commands and scope.

Signature overlap is available in `revision_signature_overlap_20261101/`:

```sh
python revision_signature_overlap_20261101/recompute_shared_gene_counts.py
```

This standard-library recount was executed during archive preparation and
matched the seven reference counts and all exact overlapping-symbol lists.
It uses frozen source symbols and the archived harmonisation map; it does not
rederive the primary signatures or reconstruct the single-cell count matrix.
See that module's README for original source-extraction dependencies and the
limitations of each published marker list. The finalized Table 11 is supplied with the journal supplementary materials;
the older internal draft has been excluded from this archive.
