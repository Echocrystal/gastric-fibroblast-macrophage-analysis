"""Export the primary GSE183904 signatures from the bundled frozen table."""

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
FROZEN_SIGNATURES = HERE / "GSE183904_PLAT_ecosystem_bulk_signatures_frozen.tsv"
PRIMARY_SIGNATURES = {
    "sender_CXCL14_POSTN_state_axis_excluded": ("sender", "PLAT-rich fibroblast"),
    "receiver_C1QC_APOE_specific_axis_excluded": ("receiver", "C1QC+APOE+ macrophage"),
}
OUTPUT_FIELDS = (
    "state_role",
    "state_label",
    "gene_symbol",
    "signature_rank",
    "equal_scoring_weight",
    "signature_id",
    "contrast_id",
    "derivation_status",
    "scoring_note",
)
SCORING_NOTE = (
    "Equal-weight mean of available-gene within-sample percentile ranks; "
    "signature_rank is selection order, not a scoring weight."
)


def main():
    with FROZEN_SIGNATURES.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream, delimiter="\t"))

    primary_rows = []
    for signature_id, (role, label) in PRIMARY_SIGNATURES.items():
        selected = [row for row in source_rows if row["signature_id"] == signature_id]
        selected.sort(key=lambda row: int(row["signature_rank"]))
        if len(selected) != 30 or [int(row["signature_rank"]) for row in selected] != list(range(1, 31)):
            raise ValueError(f"Expected exactly 30 consecutively ranked genes for {signature_id}")
        if any(row["derivation_status"] != "frozen_before_TCGA_outcome_analysis" for row in selected):
            raise ValueError(f"Unexpected derivation status for {signature_id}")
        if any(row["axis_genes_allowed"] != "FALSE" for row in selected):
            raise ValueError(f"Axis-inclusive genes found in {signature_id}")
        primary_rows.extend(
            {
                "state_role": role,
                "state_label": label,
                "gene_symbol": row["gene_symbol"],
                "signature_rank": row["signature_rank"],
                "equal_scoring_weight": 1,
                "signature_id": signature_id,
                "contrast_id": row["contrast_id"],
                "derivation_status": row["derivation_status"],
                "scoring_note": SCORING_NOTE,
            }
            for row in selected
        )

    if len({row["gene_symbol"] for row in primary_rows}) != 60:
        raise ValueError("Primary sender and receiver signatures must contain 60 distinct genes")

    output = HERE / "GSE183904_primary_sender_receiver_signatures_30plus30.tsv"
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(primary_rows)

    print(f"Exported {len(primary_rows)} primary signature genes to {output}")


if __name__ == "__main__":
    main()
