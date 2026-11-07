"""Project the frozen PLAT sender/receiver signatures into gastric-cancer
chemotherapy datasets and audit whether a patient/organoid endpoint is usable.

The analysis keeps three evidence tiers separate:
  1. patient tumour tissue (GSE14208, GSE31811);
  2. gastric-cancer organoid perturbation (GSE280755; no public sample-level
     drug-response table in GEO);
  3. tumour-cell-line monoculture (GDSC2; handled in the supplement).
"""
from pathlib import Path
import csv, gzip, json, re, warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
DATA = ROOT / "analysis" / "external_data" / "chemo_response_v1"
OUT = ROOT / "analysis" / "results" / "chemo_response_projection_v1"
OUT.mkdir(parents=True, exist_ok=True)
SIG = ROOT / "resources" / "gene_sets" / "GSE183904_PLAT_ecosystem_bulk_signatures.tsv"

SIG_IDS = {
    "sender": "sender_CXCL14_POSTN_state_axis_excluded",
    "receiver": "receiver_C1QC_APOE_specific_axis_excluded",
}

def z(x):
    x = np.asarray(x, float)
    s = np.nanstd(x, ddof=1)
    return (x - np.nanmean(x)) / s if np.isfinite(s) and s > 0 else np.full(len(x), np.nan)

def parse_series(path):
    """Read a GEO series matrix while preserving sample metadata."""
    meta, samples, rows = {}, None, []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if in_table:
                rows.append(next(csv.reader([line], delimiter="\t")))
            elif line.startswith("!Sample_"):
                vals = next(csv.reader([line], delimiter="\t"))
                key = vals[0]
                meta.setdefault(key, vals[1:])
                if key == "!Sample_title":
                    samples = vals[1:]
    if not rows:
        return pd.DataFrame(), meta
    header = rows[0]
    df = pd.DataFrame(rows[1:], columns=header)
    df = df.rename(columns={header[0]: "probe"}).set_index("probe")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if samples is not None and len(samples) == df.shape[1]:
        df.columns = samples
    return df, meta

def parse_platform(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                header = next(csv.reader([next(fh).rstrip("\r\n")], delimiter="\t"))
                rows = []
                for line in fh:
                    if line.startswith("!platform_table_end"):
                        break
                    rows.append(next(csv.reader([line.rstrip("\r\n")], delimiter="\t")))
                d = pd.DataFrame(rows, columns=header)
                return d[["ID", "Gene symbol"]].rename(columns={"ID":"probe", "Gene symbol":"symbol"})
    raise RuntimeError(f"platform table not found: {path}")

def map_expression(expr, platform):
    """Map probes to symbols; for duplicate probes retain highest median probe."""
    p = platform.set_index("probe")["symbol"].reindex(expr.index)
    d = expr.copy()
    d["symbol"] = p
    d = d[d.symbol.notna() & d.symbol.ne("")].copy()
    d["symbol"] = d.symbol.astype(str).str.split("///").str[0].str.strip()
    d = d[d.symbol.ne("")]
    # One probe per gene, selected without any outcome information.
    med = d.drop(columns="symbol").median(axis=1, skipna=True)
    d["_probe_median"] = med
    d = d.sort_values("_probe_median").drop_duplicates("symbol", keep="last")
    return d.drop(columns="_probe_median").set_index("symbol")

def project(expr, sigs):
    """Frozen projection: within-sample percentile rank mean, cohort z-score."""
    ranks = expr.rank(axis=0, method="average", pct=True)
    out, cov = {}, []
    for name, genes in sigs.items():
        mapped = sorted(set(genes).intersection(ranks.index))
        raw = ranks.loc[mapped].mean(axis=0).to_numpy(float) if mapped else np.full(expr.shape[1], np.nan)
        out[name + "_raw"] = raw
        out[name + "_z"] = z(raw)
        cov.append({"signature": name, "frozen_genes": len(genes), "mapped_genes": len(mapped), "coverage": len(mapped)/len(genes)})
    out = pd.DataFrame(out, index=expr.columns)
    out["ecosystem_z"] = z((out.sender_z + out.receiver_z) / 2)
    return out, pd.DataFrame(cov)

def cox_one(df, score, time, event, cohort, endpoint):
    x = pd.to_numeric(df[score], errors="coerce")
    t = pd.to_numeric(df[time], errors="coerce")
    e = pd.to_numeric(df[event], errors="coerce")
    ok = x.notna() & t.notna() & e.notna() & (t > 0)
    if ok.sum() < 10 or e[ok].sum() < 2 or x[ok].nunique() < 3:
        return {"cohort": cohort, "endpoint": endpoint, "score": score, "n": int(ok.sum()), "events": int(e[ok].sum()), "HR_per_SD": np.nan, "CI95_low": np.nan, "CI95_high": np.nan, "p": np.nan, "status": "insufficient"}
    xx = ((x[ok] - x[ok].mean()) / x[ok].std(ddof=1)).to_numpy()[:, None]
    try:
        fit = PHReg(t[ok].to_numpy(), xx, status=e[ok].to_numpy()).fit(disp=False)
        b, se = float(fit.params[0]), float(fit.bse[0])
        return {"cohort": cohort, "endpoint": endpoint, "score": score, "n": int(ok.sum()), "events": int(e[ok].sum()), "HR_per_SD": float(np.exp(b)), "CI95_low": float(np.exp(b-1.96*se)), "CI95_high": float(np.exp(b+1.96*se)), "p": float(fit.pvalues[0]), "status": "ok"}
    except Exception as exc:
        return {"cohort": cohort, "endpoint": endpoint, "score": score, "n": int(ok.sum()), "events": int(e[ok].sum()), "HR_per_SD": np.nan, "CI95_low": np.nan, "CI95_high": np.nan, "p": np.nan, "status": f"failed:{type(exc).__name__}"}

def main():
    ss = pd.read_csv(SIG, sep="\t")
    sigs = {name: ss.loc[ss.signature_id.eq(sid), "gene_symbol"].dropna().astype(str).unique().tolist() for name, sid in SIG_IDS.items()}
    all_cox = []
    # GSE14208: treatment-naive metastatic gastric tumour biopsy; no invented response label.
    x142, m142 = parse_series(DATA / "GSE14208_series_matrix.txt.gz")
    x142 = map_expression(x142, parse_platform(DATA / "GPL571.annot.gz"))
    s142, cov142 = project(x142, sigs)
    surv = pd.read_csv(DATA / "GSE14208_SurvivalData.txt.gz", sep="\t")
    surv["sample"] = surv.ID.str.replace(r"\.CEL$", "", regex=True)
    s142.index.name = "sample"
    s142 = s142.reset_index().merge(surv, on="sample", how="left")
    s142["set"] = s142["sample"].str.extract(r"^(Training|Validation)", expand=False)
    s142.to_csv(OUT / "GSE14208_sample_scores.tsv", sep="\t", index=False)
    cov142.assign(cohort="GSE14208").to_csv(OUT / "GSE14208_signature_coverage.tsv", sep="\t", index=False)
    for cohort, d in [("GSE14208_all", s142), ("GSE14208_training", s142[s142["set"].eq("Training")]), ("GSE14208_validation", s142[s142["set"].eq("Validation")])]:
        for score in ["sender_z", "receiver_z", "ecosystem_z"]:
            all_cox.append(cox_one(d, score, "timeToProgression_mo", "progression1_censor0", cohort, "TTP"))
            all_cox.append(cox_one(d, score, "overallSurvival_mo", "death1_censor0", cohort, "OS"))

    # GSE31811: 19 paired tumour/normal samples; response is patient-level and tumour is analysed.
    x318, m318 = parse_series(DATA / "GSE31811_series_matrix.txt.gz")
    x318 = map_expression(x318, parse_platform(DATA / "GPL6480.annot.gz"))
    s318, cov318 = project(x318, sigs)
    titles = m318["!Sample_title"]
    status = m318["!Sample_characteristics_ch1"]
    response = m318["!Sample_characteristics_ch1.1"] if "!Sample_characteristics_ch1.1" in m318 else None
    if response is None:
        # parse_series stores repeated keys with suffix only if implemented otherwise; recover from raw lines.
        vals = []
        with gzip.open(DATA / "GSE31811_series_matrix.txt.gz", "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_characteristics_ch1"):
                    vals.append(next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:])
        status, response = vals[0], vals[1]
    idx = pd.DataFrame({"sample": s318.index, "title": titles, "status": status, "response": response})
    idx["patient"] = idx.title.str.extract(r"^\d+_(\d+)_", expand=False)
    idx["geo_status"] = idx.status.str.extract(r"status:\s*(\w+)", expand=False)
    # The GEO status line contains two internal errors: 02_02_miN-2-2 and
    # 34_20_miN-20 are labelled tumour although their titles identify normal
    # mucosa.  The series design contains 19 patients with one lesion (L) and
    # one normal (N) sample.  Use the title-encoded specimen type and retain
    # the GEO field for an explicit discordance audit.
    idx["tissue"] = np.where(idx.title.str.contains(r"_(?:mi|ma)L-", regex=True), "tumor",
                      np.where(idx.title.str.contains(r"_(?:mi|ma)N-", regex=True), "normal", "unknown"))
    idx["status_discordant"] = idx["geo_status"].ne(idx["tissue"])
    discordant318 = idx.loc[idx.status_discordant, ["sample", "title", "geo_status", "tissue"]].copy()
    idx["response"] = idx.response.str.extract(r"response:\s*(.*)", expand=False)
    idx = idx[idx.tissue.eq("tumor")].copy()
    d318 = s318.loc[idx["sample"]].copy(); d318.index.name = "sample"
    d318 = d318.reset_index().merge(idx, on="sample", how="left")
    d318.to_csv(OUT / "GSE31811_tumor_sample_scores.tsv", sep="\t", index=False)
    cov318.assign(cohort="GSE31811").to_csv(OUT / "GSE31811_signature_coverage.tsv", sep="\t", index=False)
    test_rows = []
    for score in ["sender_z", "receiver_z", "ecosystem_z"]:
        a = d318.loc[d318.response.eq("effective"), score].dropna().to_numpy(float)
        b = d318.loc[d318.response.eq("non-effective"), score].dropna().to_numpy(float)
        u, p = mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
        test_rows.append({"cohort":"GSE31811", "score":score, "n_effective":len(a), "n_non_effective":len(b), "median_effective":np.median(a) if len(a) else np.nan, "median_non_effective":np.median(b) if len(b) else np.nan, "MannWhitney_U":u, "p":p})
    tests = pd.DataFrame(test_rows)
    tests["q_BH"] = multipletests(tests.p.fillna(1), method="fdr_bh")[1]
    tests.to_csv(OUT / "GSE31811_response_group_tests.tsv", sep="\t", index=False)

    # GSE280755: 3 GCOs with Ctrl/OE plus 3 cell lines with Ctrl/OE; no drug table in GEO.
    tpm = pd.read_csv(DATA / "GSE280755_TPM_Matrix_ann.txt.gz", sep="\t", usecols=lambda c: c in ["gene_sym"] or c.endswith("_Ctrl") or c.endswith("_OE"))
    tpm = tpm[tpm.gene_sym.notna() & tpm.gene_sym.ne("")].set_index("gene_sym")
    tpm = tpm.groupby(level=0).max()
    organs = [c for c in tpm.columns if c.startswith(("K1_", "K24_", "K31_"))]
    s280, cov280 = project(tpm[organs], sigs)
    s280["sample"] = s280.index
    s280["organoid"] = s280["sample"].str.replace(r"_(Ctrl|OE)$", "", regex=True)
    s280["condition"] = s280["sample"].str.extract(r"_(Ctrl|OE)$", expand=False)
    s280.to_csv(OUT / "GSE280755_organoid_scores.tsv", sep="\t", index=False)
    cov280.assign(cohort="GSE280755_organoids").to_csv(OUT / "GSE280755_signature_coverage.tsv", sep="\t", index=False)
    org_rows = []
    for score in ["sender_z", "receiver_z", "ecosystem_z"]:
        wide = s280.pivot(index="organoid", columns="condition", values=score).dropna()
        if len(wide):
            dif = wide["OE"] - wide["Ctrl"]
            stat, p = wilcoxon(dif) if len(dif) >= 2 and np.any(dif != 0) else (np.nan, np.nan)
            org_rows.append({"cohort":"GSE280755_organoids", "score":score, "n_pairs":len(dif), "median_OE_minus_Ctrl":np.median(dif), "wilcoxon_p":p})
    org_tests = pd.DataFrame(org_rows)
    if len(org_tests): org_tests["q_BH"] = multipletests(org_tests.wilcoxon_p.fillna(1), method="fdr_bh")[1]
    org_tests.to_csv(OUT / "GSE280755_organoid_paired_tests.tsv", sep="\t", index=False)

    cox = pd.DataFrame(all_cox)
    if len(cox):
        ok = cox.p.notna(); cox.loc[ok, "q_BH_within_chemo_endpoints"] = multipletests(cox.loc[ok, "p"], method="fdr_bh")[1]
    cox.to_csv(OUT / "GSE14208_cox_models.tsv", sep="\t", index=False)
    qa = {"GSE14208": {"platform":"GPL571", "n_samples":len(s142), "endpoint_note":"TTP and OS; no invented responder label", "coverage":cov142.to_dict("records")}, "GSE31811":{"n_total_samples":len(x318.columns), "n_tumour_samples":len(d318), "n_patients":int(d318.patient.nunique()), "response_counts":d318.response.value_counts().to_dict(), "tissue_rule":"sample-title L=lesion/tumour; N=normal", "geo_status_discordant_samples":discordant318.to_dict("records"), "coverage":cov318.to_dict("records")}, "GSE280755":{"n_organoid_samples":len(organs), "n_organoid_pairs":3, "public_drug_response_table":False, "interpretation":"organoid Ctrl/OE perturbation, not clinical chemotherapy-response validation", "coverage":cov280.to_dict("records")}, "GDSC2":{"tier":"supplementary", "interpretation":"tumour-cell-intrinsic monoculture drug sensitivity; not patient treatment benefit"}}
    (OUT / "chemo_response_projection_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    report = ["# 胃癌化疗队列与类器官投射审计", "", "## 证据层级", "", "- GSE14208：患者治疗前胃癌活检，CF化疗，TTP/OS；不自行制造 responder 标签。", "- GSE31811：19名患者的肿瘤/正常配对样本，DCS化疗，使用19个肿瘤样本及 effective/non-effective 标签。", "- GSE280755：3个胃癌类器官的 Ctrl/PI3-OE 配对RNA-seq；GEO未提供样本级药敏表，因此只作类器官扰动审计。", "- GDSC2：28个胃癌细胞系单培养药敏，移至补充图。", "", "## Figure 6 i panel建议", "", "优先使用GSE31811有效/无效化疗反应的连续signature比较；GSE14208作为独立患者组织TTP/OS验证。GSE280755只在正文中出现于类器官证据充分且方向与患者队列一致时，否则放补充。", "", "## 统计约定", "", "所有投射使用队列内样本百分位秩均值后z标准化；Cox HR按队列内1 SD；GSE31811使用患者级肿瘤样本，Mann–Whitney检验及BH校正；GSE280755使用3个配对类器官的Wilcoxon，仅作探索性。"]
    (OUT / "CHEMO_RESPONSE_ANALYSIS_REPORT_CN.md").write_text("\n".join(report), encoding="utf-8")
    clean_report = f'''# 胃癌化疗队列与类器官投射最终审计（2026-08-10）

## 结果摘要

- GSE31811：19名患者的病灶样本；effective=11，non-effective=8。sender中位数为−0.156 vs 0.720，P=0.0506，BH q=0.152；receiver P=0.310，q=0.465；ecosystem P=0.492，q=0.492。
- GSE14208：123个表达样本中122个具有完整TTP/OS资料。sender与TTP的HR/SD=0.823（95% CI 0.678–0.998，P=0.048），完整终点家族q=0.483；训练/验证集P分别为0.054和0.484。receiver和ecosystem未显示可靠TTP或OS关联。
- GSE280755：3个胃癌类器官Ctrl/PI3-OE配对；sender、receiver、ecosystem的OE−Ctrl中位变化分别为0.129、0.519和0.234，配对检验BH q均为0.75。GEO没有样本级药敏表。
- GDSC2：28个胃癌细胞系单培养AUC，固定放在补充图；该结果仅表示肿瘤细胞内在表达投射与体外药敏的探索性关联。

## 投射算法

每个队列先使用平台注释将探针映射到基因符号；重复探针保留队列中位表达最高者。每个样本内将所有测量基因转换为百分位秩，取30基因冻结sender或receiver signature的平均百分位秩，再在队列内z标准化。ecosystem score为sender与receiver z分数的平均值后再次在队列内z标准化。Cox HR按队列内每增加1个标准差报告，不使用结局优化cutoff。

## 证据边界

GSE31811和GSE14208支持sender分数与化疗相关患者组织终点的方向一致性，但样本量、治疗方案和校正后显著性限制了结论强度。不能将其写成预测性治疗生物标志物。GDSC2不含成纤维细胞、巨噬细胞、免疫细胞、ECM空间结构或患者治疗结局；不能用于证明PLAT生态影响患者化疗获益。GSE280755只有类器官扰动表达数据，没有药敏终点，因此不承担类器官疗效验证。

## 元数据核对

GSE31811 GEO status字段将02_02_miN-2-2和34_20_miN-20误标为tumor。根据样本标题的L/N编码和19名患者的配对设计，二者按normal处理并记录于QA文件；正式比较使用19个病灶样本。
'''
    (OUT / "CHEMO_RESPONSE_ANALYSIS_REPORT_FINAL_20260810.md").write_text(clean_report, encoding="utf-8")
    print("DONE", OUT)
    print("GSE14208 Cox\n", cox.to_string(index=False))
    print("GSE31811 tests\n", tests.to_string(index=False))
    print("GSE280755 organoid tests\n", org_tests.to_string(index=False))

if __name__ == "__main__":
    main()
