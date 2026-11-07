# Signature overlap: frozen inputs, sources and recount

This module accompanies the seven-row candidate Supplemental Table 6 panel.
The final display reports reference-set size, shared count and shared symbols.
It omits Jaccard indices, recovery fractions and the two zero-overlap rows, as
requested by the authors. Broader exploratory comparisons are preserved only
in `comparison_sets.json` / `overlap_results.json` and the original source code.
They have not been recast as positive results.

## Recompute the displayed numbers

```sh
python revision_signature_overlap_20261101/recompute_shared_gene_counts.py
```

This standard-library script uses the archived 60 primary genes, reference
symbols, measured universe and per-symbol harmonisation log. It checks all
seven counts AND exact overlapping-symbol lists against frozen expectations.
The output has exactly the six displayed columns. The original HGNC complete
download is unnecessary for this frozen-map recount. It does not rederive the
60 primary genes or regenerate the underlying single-cell atlas.

`Reference n*` is the number of distinct harmonised reference genes retained in
the GSE183904 measured universe (26,523 harmonised symbols). All 30 genes in
each primary signature are measured. Source lists before restriction and
individual symbol mappings are retained. No new statistical test is applied.

## Provenance scripts

`prepare_comparisons.py` assembled source sets; `build_review_table.py` performed
the original harmonisation and generated a wider historical review table.
Those preserved scripts use historical project paths, additional plotting
libraries, an HGNC complete-set download and/or original source workbooks.
They are provenance sources and are NOT the portable recount command above.
`fetch_hgnc.py` documents HGNC acquisition; `symbols.py` records the exact
harmonisation rule. `sources_macrophage/extract_ma_marker_sets.py` records
author-threshold extraction from Ma's original supplementary workbook, which
can be downloaded from the source links. `make_final_overlap.mjs` and
`update_shared_genes.mjs` record workbook presentation using an external
artifact-tool runtime; that runtime is not bundled or required for recounting.

Full original papers, full source supplementary workbooks, internal Word files,
Node packages and the HGNC complete download are not redistributed. Small
extracted gene-list snapshots, source citations/URLs, extraction rules and
the needed symbol mappings are included. Original absolute source paths in
metadata are historical provenance, not execution requirements for recounting.

## Interpretation of the source sets

- Li2022: the as-deposited cluster-3 positive-marker list, not an explicitly
  scored author-frozen signature. Supplement numbering and original filtering
  inconsistencies are documented, without silently redefining the set.
- Lee2025: complete 425-gene CCL2+ fibroblast and 425-gene STAT3-activated
  macrophage signatures; their counts after restriction are 410 and 402.
- Zhang2025: the complete 20-gene C1QC+ RTM row of Supplementary Table S4.
- Ma2025: author-stated adjusted-P <=0.05 and average-logFC >=0.15 thresholds;
  two unresolved numeric/date-corrupted identifiers per set were excluded
  without guessing. The recoverable valid-symbol subsets contain 626 and 805
  genes. This study reused GSE183904, so the overlap is not independent-cohort
  validation.
- Becht2016: the official MCP-counter monocytic-lineage panel. This is a
  generic lineage panel, not a macrophage-state-specific signature.

See `frozen_final_seven_sets.json` for each citation, DOI, retrieval URL,
location, extraction rule and source-specific limitation. Retain those
attributions; repository code licensing does not override third-party terms.
