# Published macrophage comparator source audit

Prepared 2026-11-01. Existing manuscript and supplementary tables were not edited.

## Verified machine-readable sets

`verified_macrophage_sets.json` contains six sets with complete retained gene
symbols, role, DOI, source path, source URL, extraction rule, and limitations.

1. Becht et al. 2016, DOI `10.1186/s13059-016-1070-5`: official MCP-counter
   **Monocytic lineage** panel, 7 genes. This is a generic lineage comparator,
   not a macrophage-specific or TAM-specific state signature.
2. Same source: **Fibroblasts** panel, 8 genes. This is a generic fibroblast
   comparator, not a CAF-specific state signature.
3. Ma et al. 2025, DOI `10.3389/fimmu.2025.1618591`, Table S2:
   **Macro_C1QC**, 626 recoverable symbols from 628 author-thresholded rows.
4. Same source: **Macro_SPP1**, 805 recoverable symbols from 807
   author-thresholded rows. These clusters are not combined.
5. Zhang et al. 2025, DOI `10.1038/s41421-025-00811-2`, Supplementary Table S4,
   **C1QC+ RTM top 20 highly variable genes**, 20 genes. This is the preferred
   fixed C1QC-state comparator. The main paper explicitly uses Table S4 genes
   for C1QC+ RTM signature scoring.
6. Same source/table: **SPP1+ C1QCmedian RTM top 20 highly variable genes**,
   20 genes. Optional additional state comparator.

## MCP-counter provenance

The complete official `Signatures/genes.txt` was downloaded from
`ebecht/MCPcounter` on GitHub. Resolved commit:
`b6eac73e91c246fcff0bb1a5c68a816cd588fc48`.
The downloaded bytes have Git blob SHA1
`a37f70367d3dc83260e966e651c5d618dbf84c7c`, matching the official Git tree
entry. `MCPcounter_verification.json` records these checks.

## Ma extraction and limitations

Full article: Europe PMC `PMC12307293/fullTextXML`.
The supplementary workbook is unmodified `Ma2025_Table2.xlsx`, obtained from
Europe PMC's `PMC12307293/supplementaryFiles` archive (`ma_supp.zip`).
Table S2 is the GSE183904 differential-expression table; Table S3 is from the
NFHGC cohort and was not substituted for Table S2.

The original Methods state: "criteria of adjusted P value <=0.05 and average
log fold change (avgLogFC) >=0.15. The top 30 genes for each cell type were
selected as representative markers". The supplied table contains many more
than 30 rows per cluster. We applied the author's stated positive-marker
thresholds to all published cluster rows, without fitting a new signature:

| Cluster | All exported rows | Rows meeting stated thresholds | Unresolved numeric identifiers | Recoverable symbols retained |
|---|---:|---:|---:|---:|
| Macro_C1QC | 716 | 628 | 2 | 626 |
| Macro_SPP1 | 1012 | 807 | 2 | 805 |

The selected source cells B17113/B17120 (C1QC) contain numeric values
45718/45911; B17440/B17732 (SPP1) contain 45717/45907. Their original Excel
number format is built-in `numFmtId=16` (date format). The source identifiers
are corrupted. No guessed reconstruction was used. These four entries are
excluded and preserved in `Ma2025_unresolved_identifiers.json` and
`Ma2025_numeric_cell_evidence.json`. Consequently, any reported overlap or
Jaccard metric describes the **recoverable-symbol subset**, not the intact
published set. An additional restriction to a measured gene universe must be
reported separately.

Ma et al. and the present signatures both use **GSE183904**. This is a
literature comparison with a shared discovery dataset, not independent
cross-dataset validation. The 626/805 sets are source-derived cluster markers,
not verified author-frozen standalone signatures. The script
`extract_ma_marker_sets.py` reproduces the extraction using the standard
library. A separate, analyst-defined top-30-by-published-logFC sensitivity is
saved in `ma_top30_sensitivity_sets.json`; it should not be described as a
verified original frozen top-30 panel.

## Zhang extraction and visual verification

Full article: Europe PMC `PMC12219098/fullTextXML`.
The publisher's complete 41-page supplement is
`Zhang2025_Supplementary_Information.pdf`. Supplementary Table S4 is on PDF page
35 (printed page 35 of 41). Both 20-gene rows were transcribed in full and
checked against the rendered page `Zhang2025_TableS4_page35.png`.
`Zhang2025_TableS4_fixed_panels.tsv` stores per-gene provenance and
`Zhang2025_TableS4_excerpt.txt` preserves the extracted table text.

An earlier temporary 5-gene RTM characteristic-marker panel was superseded
by the original scored 20-gene Table S4 panel and removed from the verified
JSON. No unspecified generic TAM gene list was invented.
