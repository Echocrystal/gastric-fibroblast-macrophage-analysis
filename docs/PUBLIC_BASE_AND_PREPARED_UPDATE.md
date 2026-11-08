# Public base and prepared supplementary resource update

Prepared locally on 8 November 2026. This folder is an update overlay for the public repository, not evidence that these additions have already been published there.

Repository: https://github.com/Echocrystal/gastric-fibroblast-macrophage-analysis

## Verified public versions

- Release `v1.0.0`: `d3be5c7a2253b4802548613626a492b83a3f8e9e`.
- Table 11 resource addition: `d19a68ba00b1fcaeb6822af1268cb41ff8c2be42`.
- Complete frozen-input clinical model reporting, and public main HEAD checked on 8 November: `c0936462c87442279844f75a0a47f3e84baca8f8`.

The later commits add files; they do not change the files in the preceding release. The checked public HEAD contains 372 files. New ECM-direction, spot-covariate and model-inclusion resources in this prepared update were absent from those public versions. They are supplied with the revised submission and remain a prepared repository update until separately published.

## Using this overlay and its manifests

Start from the public clinical reporting commit above, then copy this folder's contents to the repository root, preserving relative paths and file bytes. `UPLOAD_FILE_MANIFEST.tsv` indexes the resulting merged file tree and excludes itself to avoid a recursive self hash. It therefore includes unchanged public files that are not repeated in this smaller overlay. Its expected public base is the immutable commit above, not a future mutable main branch.

- `docs/ECM_SPATIAL_RESOURCES.md` describes fixed ECM and spatial resources, original-source provenance, requirements for separately obtained GEO matrices and bounded verification.
- `docs/ECM_SPATIAL_RESOURCE_MANIFEST.tsv` indexes the new ECM/spatial delivery.
- `resources/reporting_checks/README.md` describes complete-case denominators, original specimen/model inclusion, covariate availability and checks using frozen values.
- `scripts/verify_existing_reporting_checks.py` verifies the reporting resources without refitting models.

The original `code/CODE_FILE_MANIFEST.tsv` remains the 103-file historical source-code index; all its entries matched the checked public bytes. New source files and reproduction adapters are indexed by the new resource and merged upload manifests. Neither the historical index nor the prepared archive establishes a raw-data-to-all-figures clean-clone run.

## Byte-manifest corrections

Public `resources/RESOURCE_FILE_MANIFEST.tsv` had 10 of 11 entries recorded against original CRLF bytes; clinical `INPUT_FILE_MANIFEST.json` had 3 of 4 and `PACKAGE_FILE_MANIFEST.json` had 14 of 16. These entries have been corrected to the exact public LF byte sizes and SHA256 values, with the corrected input-manifest hash propagated into the package manifest. The historical 348-entry root upload manifest had 176 line-ending-only mismatches; the refreshed merged upload manifest uses the actual checked bytes.

All mismatches described here were verified as line-ending differences. No clinical numerical input, existing model coefficient, historical analysis script or scientific result was changed by these manifest corrections. The checks do not refit clinical or survival models.

The updated `.gitattributes` uses LF for text by default. Verbatim original-source copies and frozen ECM/spatial resources have explicit `-text` exceptions, preserving source-byte hashes, including original CRLF files. Generated files in those directories already use their recorded bytes. Avoid text-editor or transfer transformations of frozen files before checking their hashes.

The checked public README retained historical pre-upload wording. The prepared root README now corrects that wording with the verified repository/release/commit identifiers, the prepared-versus-public distinction and the clinical/new-resource entry points. This correction is local to the prepared update; no remote README edit or new publication is claimed.
