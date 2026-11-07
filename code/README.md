# WJGO computational code staging

This is a review staging copy of selected R and Python source files. It is not
yet an executable public release. Source paths are preserved relative to the
research project; `upstream_qc/` contains a byte-identical copy of the
GSE183904 quality-control script from the adjacent source project. Every code
file and its SHA-256 digest is listed in `CODE_FILE_MANIFEST.tsv`.

## Manifest status

- `current_source_identified`: one of the seven primary scripts matched to
  the current main figures. The final WJGO TIFFs also underwent journal-specific
  figure postprocessing; those scripts are staged separately.
- `candidate_current`: selected upstream analysis, figure dependency, or WJGO
  postprocessing code that appears relevant to current results. Its complete
  input/output dependency chain still needs verification before publication.
- `supporting_unverified`: an exploratory or versioned supplementary script
  retained for provenance, with final use or exact revision status unresolved.

## Release work still required

The scripts retain Windows absolute paths and references to a second local
project. Most analysis steps read processed matrices already on disk, and no
single download-and-run entry point exists. Before public release, replace
those paths with configurable project and data roots, document GEO/TCGA input
accessions and download steps, record R/Python and package versions, and check
the final figure-to-script and supplementary-figure numbering map. The staged
code contains no raw matrices, runtime environments, manuscripts, or reviewer
files. No execution result is implied by its presence here.
