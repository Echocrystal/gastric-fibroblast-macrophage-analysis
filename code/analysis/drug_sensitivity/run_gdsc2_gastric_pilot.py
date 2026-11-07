"""GDSC2 gastric-cell-line drug-sensitivity pilot for the frozen ecosystem axis.

This is deliberately a feasibility analysis.  The GDSC2 cell lines contain
tumour-intrinsic expression only; fibroblast/macrophage signatures are scored
as an exploratory transcriptomic projection and are not treated as a patient
microenvironment measurement.
"""
from pathlib import Path
import itertools
import json
import warnings

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\analysis\drug_sensitivity")
RAW = ROOT / "raw_20260807"
EXPR_DIR = RAW / "cmp_expression_by_model"
OUT = ROOT / "results_20260807"
OUT.mkdir(parents=True, exist_ok=True)

SIG_FILE = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\resources\gene_sets\GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
GDSC_FILE = RAW / "GDSC2_fitted_dose_response_27Oct23.xlsx"

GROUP_ORDER = ["Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high"]
SELECTED = {
    "5-Fluorouracil": "fluoropyrimidine",
    "Oxaliplatin": "platinum",
    "Cisplatin": "platinum",
    "Paclitaxel": "taxane",
    "Docetaxel": "taxane",
    "Irinotecan": "topoisomerase-I",
    "Gemcitabine": "pyrimidine-antimetabolite",
    "Lapatinib": "HER2/EGFR exploratory",
    "Vinorelbine": "vinca exploratory",
}

def zscore(x):
    x = np.asarray(x, float)
    s = np.nanstd(x, ddof=1)
    return (x - np.nanmean(x)) / s if np.isfinite(s) and s > 0 else np.full_like(x, np.nan)

def rank_score(expr, genes):
    genes = [g for g in genes if g in expr.index]
    if not genes:
        return np.full(expr.shape[1], np.nan), 0
    # Match the frozen bulk projection: within-sample rank / all measured genes,
    # followed by across-sample z-scoring of the signature rank mean.
    ranks = expr.rank(axis=0, method="average", pct=True)
    raw = ranks.loc[sorted(set(genes))].mean(axis=0).to_numpy(float)
    return zscore(raw), len(set(genes))

def main():
    sig = pd.read_csv(SIG_FILE, sep="\t")
    sigs = {k: g.gene_symbol.dropna().astype(str).unique().tolist()
            for k, g in sig.groupby("signature_id")}
    required = {
        "sender_primary": "sender_CXCL14_POSTN_state_axis_excluded",
        "sender_inclusive": "sender_CXCL14_POSTN_PLAT_state_inclusive",
        "receiver_specific": "receiver_C1QC_APOE_specific_axis_excluded",
    }
    files = sorted(EXPR_DIR.glob("SIDM*_rnaseq.csv"))
    if len(files) != 28:
        raise RuntimeError(f"Expected 28 expression files, found {len(files)}")
    matrices, model_meta, coverage = [], [], []
    for f in files:
        d = pd.read_csv(f, usecols=["model_id", "model_name", "gene_symbol", "rsem_tpm", "duplicate"])
        d = d[~d["duplicate"].fillna(False)].dropna(subset=["gene_symbol"])
        # A gene symbol can occur more than once; retaining the maximum TPM is
        # conservative for the small cell-line pilot.
        d = d.groupby("gene_symbol", as_index=True)["rsem_tpm"].max().to_frame()
        matrices.append(d["rsem_tpm"].rename(f.stem.replace("_rnaseq", "")))
        model_meta.append({"model_id": f.stem.replace("_rnaseq", ""), "cell_line": d.index.name})
    expr = pd.concat(matrices, axis=1).fillna(0.0)
    # Recover model names from one row per file rather than relying on gene index.
    names = []
    for f in files:
        x = pd.read_csv(f, usecols=["model_name"]).iloc[0, 0]
        names.append(str(x))
    model_ids = [f.stem.replace("_rnaseq", "") for f in files]
    sample_names = dict(zip(model_ids, names))
    score_df = pd.DataFrame({"model_id": model_ids, "cell_line": names})
    for out_name, sid in required.items():
        score, mapped = rank_score(expr, sigs[sid])
        score_df[out_name + "_z"] = score
        coverage.append({"signature": out_name, "frozen_genes": len(sigs[sid]),
                         "mapped_genes": mapped, "coverage": mapped / len(sigs[sid])})
    score_df["ecosystem_z"] = zscore((score_df.sender_primary_z + score_df.receiver_specific_z) / 2)
    sender_cut = float(score_df.sender_primary_z.median())
    receiver_cut = float(score_df.receiver_specific_z.median())
    score_df["configuration"] = np.select(
        [(score_df.sender_primary_z >= sender_cut) & (score_df.receiver_specific_z >= receiver_cut),
         (score_df.sender_primary_z >= sender_cut) & (score_df.receiver_specific_z < receiver_cut),
         (score_df.sender_primary_z < sender_cut) & (score_df.receiver_specific_z >= receiver_cut)],
        ["Dual-high", "Fibroblast-dominant", "Macrophage-dominant"], default="Dual-low")
    score_df["configuration"] = pd.Categorical(score_df.configuration, categories=GROUP_ORDER, ordered=True)
    score_df.to_csv(OUT / "GDSC2_gastric_cell_line_ecosystem_scores.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(OUT / "GDSC2_signature_coverage.tsv", sep="\t", index=False)
    score_df.groupby("configuration", observed=False).size().reindex(GROUP_ORDER).rename("n").to_csv(OUT / "GDSC2_configuration_counts.tsv", sep="\t")

    response = pd.read_excel(GDSC_FILE, usecols=["CELL_LINE_NAME", "SANGER_MODEL_ID", "CANCER_TYPE", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME", "AUC", "LN_IC50", "Z_SCORE"])
    response = response[response.CANCER_TYPE.eq("Gastric Carcinoma")].copy()
    response = response[response.DRUG_NAME.isin(SELECTED)].copy()
    response = response.merge(score_df, left_on="SANGER_MODEL_ID", right_on="model_id", how="inner")
    # Consolidate duplicate screening records for the same model and drug name.
    response = (response.groupby(["SANGER_MODEL_ID", "CELL_LINE_NAME", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME",
                                  "sender_primary_z", "sender_inclusive_z", "receiver_specific_z", "ecosystem_z", "configuration"], observed=False)
                .agg(AUC=("AUC", "median"), LN_IC50=("LN_IC50", "median"), Z_SCORE=("Z_SCORE", "median"),
                     drug_ids=("DRUG_ID", lambda x: ";".join(map(str, sorted(set(x)))))
                ).reset_index())
    response["lower_AUC_more_sensitive"] = True
    response.to_csv(OUT / "GDSC2_selected_drug_response_long.tsv", sep="\t", index=False)

    group_rows, pair_rows, cont_rows = [], [], []
    for drug, g in response.groupby("DRUG_NAME", sort=False):
        g = g.dropna(subset=["AUC", "configuration"])
        groups = [g.loc[g.configuration.eq(k), "AUC"].to_numpy(float) for k in GROUP_ORDER]
        groups_nonempty = [v for v in groups if len(v) >= 2]
        kw_stat, kw_p = (kruskal(*groups_nonempty) if len(groups_nonempty) >= 2 else (np.nan, np.nan))
        rho, sp_p = spearmanr(g["ecosystem_z"], g["AUC"]) if g.ecosystem_z.nunique() > 2 else (np.nan, np.nan)
        row = {"drug": drug, "category": SELECTED[drug], "n_models": len(g), "n_groups_ge2": len(groups_nonempty),
               "Kruskal_Wallis_H": kw_stat, "Kruskal_Wallis_p": kw_p, "Spearman_rho_ecosystem_vs_AUC": rho,
               "Spearman_p": sp_p, "median_AUC": g.AUC.median(), "median_LN_IC50": g.LN_IC50.median()}
        for k, vals in zip(GROUP_ORDER, groups):
            row[f"n_{k}"] = len(vals)
            row[f"median_AUC_{k}"] = np.nanmedian(vals) if len(vals) else np.nan
        group_rows.append(row)
        for a, b in itertools.combinations(GROUP_ORDER, 2):
            va, vb = g.loc[g.configuration.eq(a), "AUC"].to_numpy(float), g.loc[g.configuration.eq(b), "AUC"].to_numpy(float)
            if len(va) >= 2 and len(vb) >= 2:
                u, p = mannwhitneyu(va, vb, alternative="two-sided")
                pair_rows.append({"drug": drug, "group_a": a, "group_b": b, "n_a": len(va), "n_b": len(vb), "U": u, "p": p,
                                  "median_AUC_a": np.median(va), "median_AUC_b": np.median(vb)})
    groups_df = pd.DataFrame(group_rows)
    if not groups_df.empty:
        groups_df["q_Kruskal_Wallis"] = multipletests(groups_df.Kruskal_Wallis_p.fillna(1), method="fdr_bh")[1]
        groups_df["q_Spearman"] = multipletests(groups_df.Spearman_p.fillna(1), method="fdr_bh")[1]
    pairs_df = pd.DataFrame(pair_rows)
    if not pairs_df.empty:
        pairs_df["q_BH_across_drugs_pairs"] = multipletests(pairs_df.p.fillna(1), method="fdr_bh")[1]
    groups_df.to_csv(OUT / "GDSC2_selected_drug_group_tests.tsv", sep="\t", index=False)
    pairs_df.to_csv(OUT / "GDSC2_selected_drug_pairwise_tests.tsv", sep="\t", index=False)
    manifest = {"n_expression_models": len(files), "n_response_rows": len(response), "selected_drugs": SELECTED,
                "score_method": "within-cell-line percentile rank mean over frozen GSE183904 signatures; across-line z-score; ecosystem=z((sender_primary+receiver_specific)/2)",
                "interpretation": "exploratory tumour-intrinsic cell-line projection; not patient treatment response"}
    (OUT / "pilot_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("SCORES")
    print(score_df[["model_id", "cell_line", "sender_primary_z", "receiver_specific_z", "ecosystem_z", "configuration"]].to_string(index=False))
    print("\nCOVERAGE")
    print(pd.DataFrame(coverage).to_string(index=False))
    print("\nDRUG TESTS")
    print(groups_df.sort_values("Kruskal_Wallis_p").to_string(index=False))

if __name__ == "__main__":
    main()
