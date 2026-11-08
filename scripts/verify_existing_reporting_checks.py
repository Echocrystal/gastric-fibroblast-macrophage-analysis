"""Check existing-model reporting snapshots without fitting statistical models.

Run from any directory with Python 3; only the standard library is required.
The input files are resolved relative to this script's package root.
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'resources' / 'reporting_checks'


def read(name):
    with (DATA / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def truth(value):
    assert value in {'TRUE', 'FALSE', 'True', 'False'}, value
    return value.lower() == 'true'


def missing(value):
    return value in {'', 'NA'}


def rank(matrix, tol=1e-9):
    # Independent elimination on the exported designs, not an imported R rank.
    a = [list(map(float, row)) for row in matrix]
    nrows, ncols = len(a), len(a[0])
    lead = 0
    for col in range(ncols):
        pivot = max(range(lead, nrows), key=lambda i: abs(a[i][col]))
        if abs(a[pivot][col]) <= tol:
            continue
        a[lead], a[pivot] = a[pivot], a[lead]
        base = a[lead][col]
        a[lead] = [x / base for x in a[lead]]
        for i in range(lead + 1, nrows):
            multiple = a[i][col]
            a[i] = [x - multiple * y for x, y in zip(a[i], a[lead])]
        lead += 1
        if lead == nrows:
            break
    return lead


def main():
    covs = read('TCGA_CAMERA_frozen_covariates_412.tsv')
    assert len(covs) == 412
    ids = [r['sample_submitter_id'] for r in covs]
    cases = [r['case_submitter_id'] for r in covs]
    assert len(set(ids)) == len(set(cases)) == 412
    assert [int(r['sample_order']) for r in covs] == list(range(1, 413))
    by_id = {r['sample_submitter_id']: r for r in covs}
    inclusion = read('TCGA_CAMERA_patient_inclusion.tsv')
    summary = {r['model_tier']: r for r in read('TCGA_CAMERA_model_denominators.tsv')}
    assert len(inclusion) == 824
    expected = {'purity_adjusted': (400, 3), 'lineage_subtype_adjusted': (269, 8)}
    tier_ids = {}
    for tier, (n, expected_rank) in expected.items():
        rows = [r for r in inclusion if r['model_tier'] == tier]
        assert [r['sample_submitter_id'] for r in rows] == ids
        variables = summary[tier]['model_variables'].split(';')
        actual_ids = []
        for r in rows:
            c = by_id[r['sample_submitter_id']]
            expected_missing = [v for v in variables if missing(c[v])]
            assert truth(r['included']) == (not expected_missing)
            assert r['missing_model_variables'] == ';'.join(expected_missing)
            if truth(r['included']):
                actual_ids.append(r['sample_submitter_id'])
        tier_ids[tier] = actual_ids
        assert len(actual_ids) == n == int(summary[tier]['included_n'])
        assert 412 - n == int(summary[tier]['excluded_n'])
        assert n == int(summary[tier]['original_limma_reported_n'])
        design = read(f'TCGA_CAMERA_{tier}_design_matrix.tsv')
        assert [r['sample_submitter_id'] for r in design] == actual_ids
        design_cols = [k for k in design[0] if k not in {'sample_order', 'sample_submitter_id', 'case_submitter_id'}]
        assert len(design_cols) == expected_rank
        assert all(math.isfinite(float(r[k])) for r in design for k in design_cols)
        assert rank([[r[k] for k in design_cols] for r in design]) == expected_rank
        assert n - expected_rank == int(summary[tier]['residual_df'])
        for r in design:
            c = by_id[r['sample_submitter_id']]
            assert float(r['(Intercept)']) == 1
            for k in design_cols:
                if k in c:
                    assert math.isclose(float(r[k]), float(c[k]), rel_tol=1e-12, abs_tol=1e-12)
                elif k.startswith('tcga2014_subtype'):
                    assert float(r[k]) == int(c['tcga2014_subtype'] == k.removeprefix('tcga2014_subtype'))
    assert set(tier_ids['lineage_subtype_adjusted']) <= set(tier_ids['purity_adjusted'])
    patterns = Counter((missing(r['purity_z_model']), missing(r['tcga2014_subtype'])) for r in covs)
    assert patterns == {(False, False): 269, (False, True): 131, (True, False): 5, (True, True): 7}
    for name in ['ecosystem_z', 'fibroblast_z_model', 'monocytic_z_model']:
        assert all(not missing(r[name]) for r in covs)
    for name in ['purity_z_model', 'fibroblast_z_model', 'monocytic_z_model']:
        values = [float(r[name]) for r in covs if not missing(r[name])]
        assert abs(statistics.mean(values)) < 1e-12
        assert abs(statistics.stdev(values) - 1) < 1e-12
    hash_checks = read('TCGA_CAMERA_original_input_hash_checks.tsv')
    assert len(hash_checks) == 3
    assert all(r['current_sha256'] == r['original_analysis_sha256'] and truth(r['exact_match']) for r in hash_checks)
    # Cross-check each survival inclusion list against its original-result n.
    survival_counts = defaultdict(Counter)
    for r in read('TCGA_survival_patient_inclusion.tsv'):
        key = (r['analysis_family'], r['exposure'], r['tier'])
        survival_counts[key][truth(r['included'])] += 1
        assert truth(r['included']) == (r['exclusion_reason'] == '')
    survival = read('TCGA_survival_denominator_checks.tsv')
    assert len(survival) == len(survival_counts) == 21
    for r in survival:
        key = (r['analysis_family'], r['exposure'], r['tier'])
        assert survival_counts[key][True] == int(r['included_n']) == int(r['original_reported_n'])
        assert sum(survival_counts[key].values()) == int(r['candidate_n']) == 412
        assert int(r['events']) == int(r['original_reported_events'])
    kim = read('Kim2018_analysis_input_field_availability.tsv')
    kim_fields = {r['variable']: r for r in kim}
    for name in ['treatment_line', 'tumour_anatomic_site', 'histology']:
        r = kim_fields[name]
        assert not truth(r['source_field_present'])
        assert int(r['available_n']) == 0 and int(r['unavailable_n']) == 45
        assert r['availability_type'] == 'field_not_provided_in_analysed_clinical_input'
    assert int(kim_fields['msi_type_pentaplex']['available_n']) == 45
    assert int(kim_fields['ebv_in_situ']['available_n']) == 45
    assert int(kim_fields['pdl1_cps']['available_n']) == 41
    assert int(kim_fields['pdl1_cps']['unavailable_n']) == 4
    fits = read('Kim2018_stored_fit_denominator_checks.tsv')
    assert len(fits) == 11
    assert all(int(r['stored_fit_n']) == 45 and int(r['responders']) == 12 and int(r['missing_or_excluded_n']) == 0 for r in fits)
    kim_inclusion = read('Kim2018_existing_model_patient_inclusion.tsv')
    assert len(kim_inclusion) == 495 and all(truth(r['included']) for r in kim_inclusion)
    assert len(set(r['patient_id'] for r in kim_inclusion)) == 45
    assert all(r['patient_id'].startswith('KIM') for r in kim_inclusion)
    return {
        'pass': True, 'scope': 'Internal reporting-resource checks; no model fitting.',
        'TCGA_CAMERA_n': {tier: len(v) for tier, v in tier_ids.items()},
        'TCGA_CAMERA_design_ranks': {tier: v[1] for tier, v in expected.items()},
        'TCGA_survival_original_n_event_checks': len(survival),
        'Kim2018_stored_fit_denominator_checks': len(fits),
        'Kim2018_patients_per_model': 45, 'Kim2018_responders_per_model': 12,
        'Kim2018_absent_fields_verified': ['treatment_line', 'tumour_anatomic_site', 'histology'],
    }


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
