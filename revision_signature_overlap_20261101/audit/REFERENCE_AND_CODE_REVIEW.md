# Independent review of reference sets and calculation code

Reviewed 2026-11-01: `build_review_table.py`, `symbols.py`, the specific CAF set JSON files, and `verified_macrophage_sets.json`. No manuscript or root build script was changed by this review.

The overlap calculation is correct for the explicitly defined, harmonized GSE183904 feature universe. Both frozen primary signatures retain 30/30 genes. The original 26,571 unique exported symbols map to 26,523 unique normalized identifiers. The code properly deduplicates each reference after harmonization, calculates the measured-space intersection and union, and uses recovery k/30. The code labels signatures, marker panels, ECM panels and source-derived positive marker sets separately. It explicitly acknowledges that Ma et al. share GSE183904 and therefore are not independent validation.

## Material issues found

1. **Date-corrupted Ma identifiers.** The original extracted C1QC list contained numeric identifiers `45718` and `45911`; SPP1 contained `45717` and `45907`. The source agent checked the original XLSX and confirmed numeric/date cell formatting (numFmt16). The extraction agent removed these unresolved identifiers, retained a log, and did not guess gene names. The updated C1QC set has 626 recoverable symbols from 628 selected raw entries; SPP1 has 805 recoverable symbols from 807 selected raw entries. Final labels and interpretation must say recoverable-symbol subset. These are not intact original sets. Excluding them does not change the known overlap because all primary receiver genes are valid and already recovered.

2. **Column label ambiguity.** `original_n = len(B)` is the harmonized and deduplicated full-reference size. It is not necessarily the raw published entry count. Prefer `Full reference set n (harmonised)`; if a published/raw count is shown, use `raw_unique_n` and separately report unresolved symbols.

3. **PDF width check passed.** The inspected PDF column widths total **1,109 pt**, within the **1,134.55 pt** usable landscape-A3 width after 28-pt left and right margins. An earlier reviewer message incorrectly added the widths as 1,209 pt; that issue is withdrawn. Final visual inspection remains appropriate.

## Independently recalculated current reference metrics

| Reference | Valid full set n | In M | Overlap k | Jaccard | Recovery |
|---|---:|---:|---:|---:|---:|
| Lee CCL2 fibroblast | 425 | 410 | 13 | 0.030444965 | 13/30 |
| Lee STAT3 macrophage | 425 | 402 | 9 | 0.021276596 | 9/30 |
| Li eCAF cluster 3 as deposited | 263 | 263 | 5 | 0.017361111 | 5/30 |
| Li eCAF significance-only sensitivity | 261 | 261 | 5 | 0.017482517 | 5/30 |
| Yang ECM hub panel | 10 | 10 | 0 | 0 | 0/30 |
| MCP-counter monocytic lineage | 7 | 7 | 1 | 0.027777778 | 1/30 |
| MCP-counter fibroblasts | 8 | 8 | 0 | 0 | 0/30 |
| Zhang C1QC RTM fixed top20 panel (final selection) | 20 | 20 | 1 | 0.020408163 | 1/30 |
| Ma C1QC recoverable-marker subset | 626 | 626 | 30 | 0.047923323 | 30/30 |
| Ma SPP1 recoverable-marker subset | 805 | 805 | 17 | 0.020782396 | 17/30 |

The initial source JSONs were independently recalculated by `audit_reference_sets.py`, with results in `reference_set_independent_audit.json`. The **final** 11-row output (9 computable rows) was independently checked by `audit_final_output.py`; see `FINAL_OUTPUT_AUDIT.json`. That final check passed all numerical calculations, Excel cell values, existence of all 9 source files, 11 source hashes, and unchanged hashes for all 5 original submission inputs. Zhang's final comparator is the author-published Table S4 top20 panel, replacing the initial five characteristic RTM markers. The Li significance-only sensitivity list is not an additional main-table row.

One final source-location wording correction was requested: Zhang's `source_location` retained the older `Methods: Definitions of RTMs` label, while its correct extraction rule and source file already point to Supplementary Table S4, PDF page 35. The final table should update this location label.

## Interpretation

The receiver set is fully contained in the recoverable Ma C1QC positive-marker list. This directly supports an overlap statement, but it compares signatures drawn from the same discovery dataset and a broad marker list; it does not establish independent reproducibility or identity of the scoring signatures. The low Jaccard values alongside high recovery reflect reference set size. Neither zero overlap with a five-gene marker panel nor zero overlap with a generic lineage panel proves that a cell state is novel.

Suggested concise panel title: **Signature overlap with published CAF and macrophage gene sets**. The set class column and footnotes can retain the exact distinction among signatures, markers and ECM programmes.
