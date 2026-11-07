> Historical assessment of the original staged code. Later revision modules
> and current execution scope are described in CURRENT_REVISION_MAP.md.

# Release readiness

## Complete in local staging

- Current main-figure script mapping and copies of selected analysis scripts.
- Supplemental Figure 1-13 mapping; 11 script sources identified, one
  candidate (Figure 5), and one unverified final variant (Figure 11).
- Frozen primary signature list and state marker source tables.
- Eleven additional frozen gene/program and ligand-receptor TSV resources,
  totaling 1,289,687 bytes, with source paths, hashes, and attribution notes
  in `resources/RESOURCE_FILE_MANIFEST.tsv` and `resources/README.md`.
- Seed audit with script locations and manuscript-method coverage.
- Public input accession list and file hashes for copied code.

## Static checks completed

- All 103 staged code files matched their entries in
  `code/CODE_FILE_MANIFEST.tsv` by SHA-256.
- All 11 staged resource TSVs matched both their original project files and
  `resources/RESOURCE_FILE_MANIFEST.tsv` by byte count and SHA-256.
- Python AST parsing succeeded for all 47 staged `.py` files, and R `parse()`
  succeeded for all 56 staged `.R` files. These checks only establish file
  integrity and syntax; they did not execute any analysis or verify outputs.

## Work required before claiming a portable rerun

- Path portability: multiple scripts contain hard-coded `I:/` or `I:\\`
  locations, and some refer to a separate historical source project.
- Input acquisition: raw public accessions are named, but many scripts read
  already processed matrices and there is no unified download entrypoint.
- Environment capture: the manuscript names the principal Python and R
  versions, but no complete dependency lockfile has been validated against
  the selected scripts.
- End-to-end check: selected scripts have not been rerun in a clean clone of
  this staging directory. Figure typography was also adjusted for WJGO after
  analytic figure generation.
- Resource acquisition: the original NicheNet-v2 RDS and processed matrices
  are not included in `resources/`; source and reuse terms for the frozen
  derived TSVs are documented in `resources/README.md`.
- Editorial integration: the manuscript's Data sharing statement does not
  yet contain a code URL. A public URL must not be inserted before release.

This is a code-availability package under preparation. It is already suitable
for internal review of script provenance, but the reviewer should not be told
that a fresh clone reproduces every figure until the remaining work is done.
