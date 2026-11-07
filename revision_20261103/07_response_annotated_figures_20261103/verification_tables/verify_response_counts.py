"""Independently aggregate the frozen 45-patient table; no hypothesis tests."""
from pathlib import Path
import csv
import hashlib
import json
import math
from collections import Counter
from statistics import NormalDist

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
SOURCE = ROOT / '02_clinical_stratification/Kim2018_patients_with_restriction_flags.tsv'
CONFIGS = ['Dual-low', 'Fibroblast-dominant', 'Macrophage-dominant', 'Dual-high']
Z = NormalDist().inv_cdf(0.975)


def write(name, rows):
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def wilson(k, n):
    if not n:
        return '', ''
    p = k / n
    den = 1 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / den
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return max(0.0, center - half), min(1.0, center + half)


def summary(rows):
    n = len(rows)
    k = sum(int(x['response_binary']) for x in rows)
    counts = Counter(x['best_response'] for x in rows)
    return dict(n=n, responders_n=k, nonresponders_n=n-k,
                response_fraction=k/n if n else '',
                CR_n=counts['CR'], PR_n=counts['PR'],
                SD_n=counts['SD'], PD_n=counts['PD'])


with SOURCE.open(encoding='utf-8-sig', newline='') as f:
    patients = list(csv.DictReader(f, delimiter='\t'))
assert len(patients) == 45
assert len({x['patient_id'] for x in patients}) == 45
assert set(x['ecotype4'] for x in patients) == set(CONFIGS)
expected = {'CR': 1, 'PR': 1, 'SD': 0, 'PD': 0}
for x in patients:
    assert expected[x['best_response']] == int(x['response_binary']), x['patient_id']
    assert x['msi_type_pentaplex'] in ['MSI', 'MSS']
    assert x['ebv_in_situ'] in ['Positive', 'Negative']
    assert int(x['mss_and_ebv_negative_flag']) == int(
        x['msi_type_pentaplex'] == 'MSS' and x['ebv_in_situ'] == 'Negative')
    cps_is_missing = x['pdl1_cps'].strip().lower() in ['', 'nan', 'na', 'none']
    assert int(x['cps_available_flag']) == int(not cps_is_missing)

for col, levels, filename in [
    ('msi_type_pentaplex', ['MSI', 'MSS'], 'configuration_MSI_by_response.tsv'),
    ('ebv_in_situ', ['Positive', 'Negative'], 'configuration_EBV_by_response.tsv'),
]:
    out = []
    for config in CONFIGS:
        for level in levels:
            subset = [x for x in patients if x['ecotype4'] == config and x[col] == level]
            out.append(dict(configuration=config, biomarker=col, status=level, **summary(subset)))
    assert sum(x['n'] for x in out) == 45
    assert sum(x['responders_n'] for x in out) == 12
    write(filename, out)

response_rows = []
for population in ['All patients', 'MSS and EBV ISH-negative']:
    for config in CONFIGS:
        subset = [x for x in patients if x['ecotype4'] == config and (
            population == 'All patients' or
            (x['msi_type_pentaplex'] == 'MSS' and x['ebv_in_situ'] == 'Negative'))]
        s = summary(subset)
        low, high = wilson(s['responders_n'], s['n'])
        response_rows.append(dict(population=population, configuration=config, **s,
                                  Wilson95_low=low, Wilson95_high=high,
                                  CI_method='Two-sided Wilson score interval; z=Phi^-1(0.975)'))
write('full_and_MSS_EBVnegative_response_Wilson95.tsv', response_rows)

missing = [x for x in patients if not int(x['cps_available_flag'])]
assert len(missing) == 4
write('CPS_missing_patients.tsv', [dict(
    patient_id=x['patient_id'], configuration=x['ecotype4'],
    MSI_status=x['msi_type_pentaplex'], EBV_ISH_status=x['ebv_in_situ'],
    best_response=x['best_response'], response_binary=x['response_binary'],
    CPS_missing=True) for x in missing])
write('CPS_missing_response_summary.tsv', [
    dict(configuration=config, **summary(missing if config == 'All configurations' else
        [x for x in missing if x['ecotype4'] == config]))
    for config in CONFIGS + ['All configurations']
])
write('best_response_binary_validation.tsv', [dict(
    best_response=label, expected_response_binary=expected[label],
    observed_n=sum(x['best_response'] == label for x in patients),
    observed_binary_1_n=sum(x['best_response'] == label and x['response_binary'] == '1' for x in patients),
    mismatches_n=sum(x['best_response'] == label and int(x['response_binary']) != expected[label] for x in patients)
) for label in ['CR', 'PR', 'SD', 'PD']])

outputs = sorted(OUT.glob('*.tsv'))
qa = dict(source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          independent_patients_n=len(patients), responders_n=12,
          restricted_patients_n=sum(x['n'] for x in response_rows if x['population'] != 'All patients'),
          CPS_missing_n=len(missing), CPS_missing_responders_n=summary(missing)['responders_n'],
          response_mapping='CR/PR=1; SD/PD=0; all 45 verified',
          hypothesis_tests_performed=False, new_models_fitted=False,
          output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs})
(OUT / 'verification_QA.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(qa, ensure_ascii=False, indent=2))
