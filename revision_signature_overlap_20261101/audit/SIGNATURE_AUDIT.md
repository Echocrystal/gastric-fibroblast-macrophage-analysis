# Signature and overlap audit — 2026-11-01

Primary signature source: `I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\WJGOsubmit\11_Code_Release_Staging\gene_sets\GSE183904_primary_sender_receiver_signatures_30plus30.tsv`.

The five frozen/list source files contain exactly the same two primary signatures in the same selection order. Each primary set has 30 unique genes; neither contains PLAT, ANXA2, ITGAM, ITGB2, LRP1 or ILK. The 120-row frozen files additionally contain a PLAT-inclusive sender and a merged C1QC receiver sensitivity signature; those must not be substituted for the primary sets.

The GSE183904 feature universe contains **26,571 unique gene symbols**; both primary signatures have **30/30** represented. Source: `I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\analysis\processed\GSE183904_global_atlas\GSE183904_global_atlas_genes.tsv`, column `gene_symbol`. The export script (45, lines 71–77) obtains these rows from the raw count assay and asserts identical gene order across the 26 primary tumour samples. This is a measured feature universe, without imposing a new expression threshold.

For this explicitly defined M, Jaccard = k/(30+nB−k) and recovery = k/30 are valid. For any other M, coverage must be rechecked before using 30 in the Jaccard denominator.

The WJGO and MedComm workbooks contain the same primary 60 genes in the S6 gene annotation panel. Their original S6 panels cover clinicopathology, biomarkers, TIDE, CAMERA, MCP-counter and drug sensitivity. No dedicated published CAF/macrophage signature overlap panel is present. The existing MCP-counter overlap-removal panel records **Monocytic lineage: PLA2G7, n=1**, which is not a CAF/TAM novelty analysis.

## Exact selection rule found in script 44

Sender: fib_c1 versus fib_c4. Receiver: my_c3 versus my_c145. Patient-blocked edgeR quasi-likelihood estimates are filtered by FDR ≤0.01, logFC ≥0.5, target patient detection ≥0.5, paired direction consistency ≥0.65, TCGA gene availability, and technical exclusions. The six focal-axis genes are removed. Eligible genes are ordered by ascending FDR, descending consistency and descending logFC, then the first 30 are frozen. TCGA gene availability is used; outcome-based selection is not used in the code.

## Corrections to the supplied template

- If A is restricted to M, Jaccard denominator is |A intersect M| + |B intersect M| - overlap; the simplified 30 + n - overlap denominator is valid only if all 30 primary signature genes are in M.
- Recovery overlap/30 is an explicitly fixed original-signature denominator, not coverage-conditioned recovery. If any signature genes are absent in M, report signature coverage and optionally overlap/|A intersect M| separately.
- A generic review framework is not a fixed signature. A non-retrieved list means not computable from retrieved evidence, not proof that no gene list exists.
- Do not assert substantial overlap before computing it. A small Jaccard can coexist with appreciable recovery when the published set is large; both metrics should be interpreted with set sizes.
- The actual frozen signatures exclude all six focal-axis genes but do not necessarily contain their naming markers (CXCL14, POSTN, C1QC, APOE). Comparing only naming markers would test a different set.
- Full harmonized gene sets without expression-space restriction can be a primary set-definition comparison; measured-space restriction may be a separate sensitivity analysis if its M is precisely documented.
