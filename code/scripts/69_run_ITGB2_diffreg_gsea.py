#!/usr/bin/env python3
"""Ranked-gene-set enrichment after ITGB2 receiver-side virtual perturbation.

The analysis is deliberately limited to the six frozen Figure 7 programmes.
ITGB2 is removed from the ranking so that the self-hit cannot create artificial
programme enrichment. Each algorithmic seed is analysed independently, then
summarised by median NES and seed-direction concordance.
"""

from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
INPUT = PROJECT / "figures/source_data/Figure7_receiver_screen_v3/Figure7B_ITGB2_diffreg_all_seeds.tsv"
PROGRAMS = PROJECT / "resources/gene_sets/GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv"
OUTDIR = PROJECT / "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2/ITGB2_GSEA_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

N_PERM = 2000
MIN_SIZE = 10
MAX_SIZE = 200
BASE_SEED = 2026081101


def bh(values):
    values = np.asarray(values, dtype=float)
    n = len(values)
    order = np.argsort(values)
    ranked = values[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0, 1)
    return out


def enrichment_score(scores, gene_set):
    """Weighted signed GSEA enrichment score and leading-edge genes."""
    genes = np.asarray(scores.index)
    vals = scores.to_numpy(dtype=float)
    hit = np.isin(genes, list(gene_set))
    nh = int(hit.sum())
    if nh == 0 or nh == len(genes):
        return np.nan, []
    weights = np.abs(vals) ** 1.0
    hit_weight = weights * hit
    hit_weight_sum = hit_weight.sum()
    if hit_weight_sum <= 0:
        hit_weight = hit.astype(float)
        hit_weight_sum = hit_weight.sum()
    running = np.cumsum(hit_weight / hit_weight_sum) - np.cumsum((~hit) / (len(genes) - nh))
    imax = int(np.argmax(running))
    imin = int(np.argmin(running))
    if abs(running[imax]) >= abs(running[imin]):
        es = float(running[imax])
        leading = genes[: imax + 1][hit[: imax + 1]].tolist()
    else:
        es = float(running[imin])
        leading = genes[imin:][hit[imin:]].tolist()
    return es, leading


def permutation_null(scores, set_size, rng):
    genes = np.asarray(scores.index)
    vals = scores.to_numpy(dtype=float)
    n = len(genes)
    out = np.empty(N_PERM, dtype=float)
    # The null gene set is sampled from the frozen ranking universe, preserving
    # the observed set size while removing pathway membership.
    for i in range(N_PERM):
        selected = rng.choice(n, size=set_size, replace=False)
        hit = np.zeros(n, dtype=bool)
        hit[selected] = True
        weights = np.abs(vals) * hit
        denom = weights.sum()
        if denom <= 0:
            weights = hit.astype(float)
            denom = weights.sum()
        running = np.cumsum(weights / denom) - np.cumsum((~hit) / (n - set_size))
        out[i] = running[np.argmax(np.abs(running))]
    return out


def main():
    diff = pd.read_csv(INPUT, sep="\t")
    programs = pd.read_csv(PROGRAMS, sep="\t")
    diff["gene"] = diff["gene"].str.upper()
    programs["gene_symbol"] = programs["gene_symbol"].str.upper()
    frozen = {
        p: set(g["gene_symbol"].dropna())
        for p, g in programs.groupby("program")
    }
    rows = []
    leading_rows = []
    for seed in sorted(diff.seed.unique()):
        d = diff[diff.seed == seed].drop_duplicates("gene").copy()
        d = d[d.gene != "ITGB2"].copy()
        d = d[np.isfinite(d.Z)].sort_values("Z", ascending=False)
        scores = pd.Series(d.Z.to_numpy(float), index=d.gene.to_numpy())
        universe = set(scores.index)
        for j, (program, genes) in enumerate(frozen.items()):
            gs = genes & universe
            if not (MIN_SIZE <= len(gs) <= MAX_SIZE):
                continue
            es, leading = enrichment_score(scores, gs)
            rng = np.random.default_rng(BASE_SEED + int(seed) + j)
            null = permutation_null(scores, len(gs), rng)
            p = (1 + np.sum(np.abs(null) >= abs(es))) / (N_PERM + 1)
            nes = es / np.mean(np.abs(null)) if np.mean(np.abs(null)) > 0 else np.nan
            direction = "positive" if es > 0 else "negative"
            rows.append({
                "seed": int(seed), "program": program,
                "n_ranked_genes": len(scores), "set_size_in_universe": len(gs),
                "ES": es, "NES": nes, "permutation_p": p,
                "direction": direction, "leading_edge_n": len(leading),
            })
            for rank, gene in enumerate(leading, 1):
                leading_rows.append({"seed": int(seed), "program": program,
                                     "leading_edge_rank": rank, "gene": gene})
    result = pd.DataFrame(rows)
    result["FDR_within_seed"] = result.groupby("seed")["permutation_p"].transform(lambda x: bh(x.to_numpy()))
    result["FDR_global"] = bh(result.permutation_p.to_numpy())
    result.to_csv(OUTDIR / "ITGB2_diffreg_GSEA_by_seed.tsv", sep="\t", index=False)
    pd.DataFrame(leading_rows).to_csv(OUTDIR / "ITGB2_diffreg_GSEA_leading_edge.tsv", sep="\t", index=False)

    summary_rows = []
    for program, g in result.groupby("program"):
        summary_rows.append({
            "program": program,
            "n_seeds": len(g),
            "median_ES": g.ES.median(),
            "median_NES": g.NES.median(),
            "min_NES": g.NES.min(), "max_NES": g.NES.max(),
            "positive_seed_count": int((g.ES > 0).sum()),
            "negative_seed_count": int((g.ES < 0).sum()),
            "seed_FDR_lt_0_05_count": int((g.FDR_within_seed < 0.05).sum()),
            "median_permutation_p": g.permutation_p.median(),
            "median_leading_edge_n": g.leading_edge_n.median(),
        })
    summary = pd.DataFrame(summary_rows).sort_values("median_NES", ascending=False)
    summary.to_csv(OUTDIR / "ITGB2_diffreg_GSEA_summary.tsv", sep="\t", index=False)
    manifest = {
        "generated_at": "2026-08-11",
        "input": str(INPUT),
        "gene_set_source": str(PROGRAMS),
        "gene_set_policy": "six frozen Figure 7 programmes; no post-hoc gene-set addition",
        "ranking": "signed scTenifoldKnk differential-regulation Z-score; ITGB2 self-hit removed",
        "permutation_count": N_PERM,
        "minimum_set_size": MIN_SIZE,
        "maximum_set_size": MAX_SIZE,
        "seeds": sorted(int(x) for x in diff.seed.unique()),
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "gene_set_sha256": hashlib.sha256(PROGRAMS.read_bytes()).hexdigest(),
        "interpretation": "ranked functional-programme enrichment after virtual perturbation; not expression-level experimental GSEA",
    }
    (OUTDIR / "ITGB2_diffreg_GSEA_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(result.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
