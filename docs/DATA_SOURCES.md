# Data sources and input boundary

The study reused public datasets. This staging package does not redistribute
their raw matrices or clinical records. Retrieve each input from its source
and retain the accession or release version used.

| Role | Source | Identifier or location |
| --- | --- | --- |
| Broad single-cell discovery | GEO | [GSE206785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206785) |
| Gastric-cancer cell-state definition and signatures | GEO | [GSE183904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183904) |
| Spatial validation, main and repeated sections | GEO | [GSE270678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270678) |
| Independent spatial validation | GEO | [GSE251950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251950) |
| Bulk gastric-cancer projection | GDC/TCGA | [TCGA-STAD](https://portal.gdc.cancer.gov/projects/TCGA-STAD) |
| External bulk cohorts | GEO | [GSE62254](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62254), [GSE84437](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437), [GSE15459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15459) |
| Observed anti-PD-1 response | Kim et al., Nature Medicine 2018 | [DOI 10.1038/s41591-018-0101-z](https://doi.org/10.1038/s41591-018-0101-z) |
| Gastric-cancer cell-line expression | [Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/) | RSEM-TPM release used in the study |
| Cell-line drug response | [GDSC2](https://www.cancerrxgene.org/downloads/bulk_download) | fitted dose-response AUC file dated 27 Oct 2023 |

The original scripts generally expect normalized or otherwise processed
tables at the directory names used in the local project. Those derived input
paths are visible in the scripts. A portable download-and-preparation workflow
has not yet been validated for all modules. The manuscript's Supplementary
Methods describe source-specific filtering and the independent inferential
unit for each analysis.

The GSE183904 secondary-QC script originated in a separate upstream project.
A byte-identical source copy is staged at
`code/upstream_qc/14_prepare_GSE183904_manifest_and_secondary_qc.R`; its
ordered sample manifest is not included. That manifest and the processed input
chain must be documented before claiming a self-contained clean rerun of the
single-cell atlas.

The 60 fixed primary signature genes are in
`gene_sets/GSE183904_primary_sender_receiver_signatures_30plus30.tsv`.
The source table also includes alternate signatures; these should not be
substituted for the two primary signatures.
