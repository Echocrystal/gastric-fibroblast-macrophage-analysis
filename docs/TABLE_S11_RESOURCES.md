# Table S11 resource files

This archive accompanies Table S11. docs/TABLE_S11_RESOURCE_FILE_MANIFEST.tsv lists the files,
their contents and the related resource or analysis identifiers. All file paths
are relative to the archive root and can be retained in the study code repository.

## Contents

- resources/gene_sets/ and resources/ligand_receptor/: eleven fixed analysis
  resource TSVs, including annotation modules, programme gene sets, gene universes,
  matched controls and the ligand-receptor prior.
- gene_sets/: the primary sender and receiver signatures (30 genes each), plus
  descriptive markers for fibroblast and myeloid clusters. The full signature
  source also contains two alternative 30-gene lists. The four lists remain separate.
- resources/signature_references/: reference-list rows used in the seven Table S6
  comparisons, with original and harmonised gene symbols and membership flags.
- resources/patient_summaries/: fixed patient-state features, fractions and counts
  used by the ecology analyses, including the ECM-direction values read by the scripts.
- revision_20261103/: the Kim2018 patient input and spatial endpoint estimates,
  coverage and calibration used by the revision analyses.

Short marker and control lists are printed in S11C. Longer programme memberships
are listed gene by gene in S11D. Original numerical eligibility and matching fields
remain in the TSV files. Analysis scripts are in the separate study code package.

## Source attribution

The supplied Hallmark, GO-BP and Reactome programme subsets were prepared with
msigdbr 26.1.0. Cite MSigDB (DOI:10.1016/j.cels.2015.12.004). Source terms:
https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp .

The NicheNet-v2 human network is by Robin Browaeys, Zenodo record 7074291:
https://doi.org/10.5281/zenodo.7074291 (CC BY 4.0). Project processing uppercased
symbols, consolidated duplicate ligand-receptor pairs and source annotations,
and added benchmark and study-data availability flags. Original source annotations
remain in the TSVs. The source NicheNet RDS is available from the Zenodo record.

Kim2018 clinical input derives from the published study
(DOI:10.1038/s41591-018-0101-z). Other published reference-set DOIs are listed in
S11A and in the Table S6 references worksheet. The remaining signatures,
annotation modules, gene universes and matched controls were derived for this study.
