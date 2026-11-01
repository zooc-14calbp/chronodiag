"""Reproduce the 2026-11-02 analytical comparisons and point-only Figure 2.

Run from any directory with Python, NumPy and Matplotlib installed:
  python revision_support/reproduce_revision.py --output-dir revision_support/reproduced

Inputs are resolved within the repository or an extracted public source snapshot. No new resampling,
surrogate simulation, detrending, network request or source-data edit is made.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def close(a, b):
    if not math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-11):
        raise AssertionError((a, b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--skip-figure', action='store_true')
    args = parser.parse_args()
    root = args.repo_root.resolve()
    out = args.output_dir.resolve() if args.output_dir else root / 'revision_support/reproduced'
    if out == root or out == root / 'revision_support/verified_results':
        raise ValueError('Choose a separate reproduction output directory.')
    out.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def rows(path):
        inputs[path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
            return list(csv.DictReader(stream))

    def single(rs, dataset, **conditions):
        matches = [r for r in rs if r['dataset'] == dataset
                   and all(r[k] == str(v) for k, v in conditions.items())]
        if len(matches) != 1:
            raise AssertionError((dataset, conditions, len(matches)))
        return matches[0]

    answers, deff_rows, figure_rows = [], [], []
    for dataset, stage in [('ca533', 'itrdb_empirical_stage_b_v01'),
                           ('fran003', 'itrdb_empirical_stage_b_v01'),
                           ('yoshino', 'yoshino_empirical_stage_c_v01')]:
        base = root / 'analysis' / stage
        s = single(rows(base / 'results/structure_summary_v01.csv'), dataset,
                   detrending_method='AgeDepSpline')
        N, T, rho = int(s['n_trees']), int(s['analysis_years']), float(s['rbar'])
        eps = lambda k: k * rho / (1 + (k - 1) * rho)
        close(eps(N), s['eps'])
        n, m = N // 2, N - N // 2
        subset = single(rows(base / 'results/subset_summary_v01.csv'), dataset, subset_size=n)
        split = single(rows(base / 'results/split_half_summary_v01.csv'), dataset)
        boot = single(rows(base / 'results/bootstrap_summary_v01.csv'), dataset)
        raw = [r for r in rows(base / 'preprocessed/tree_level_all_methods_v01.csv')
               if r['dataset'] == dataset and r['detrending_method'] == 'AgeDepSpline']
        years, ids = sorted({int(r['year']) for r in raw}), sorted({r['tree_id'] for r in raw})
        if len(years) != T or len(ids) != N or len(raw) != N * T:
            raise AssertionError('Complete-cohort dimensions do not match the archived summary.')
        yi, ii = {v: i for i, v in enumerate(years)}, {v: i for i, v in enumerate(ids)}
        x = np.full((T, N), np.nan)
        for r in raw:
            pos = yi[int(r['year'])], ii[r['tree_id']]
            if np.isfinite(x[pos]):
                raise AssertionError('Duplicate year/individual observation.')
            x[pos] = float(r['value_z'])
        assert np.isfinite(x).all()
        assert np.max(np.abs(x.mean(axis=0))) < 1e-12
        assert np.max(np.abs(x.var(axis=0, ddof=1) - 1)) < 1e-12
        cov = np.cov(x, rowvar=False, ddof=1)
        full_var = float(np.var(x.mean(axis=1), ddof=1))
        close(full_var, s['chronology_variance'])
        close((cov.sum() - np.trace(cov)) / (N * (N - 1)), rho)
        expected_boot = 1 + (N - 1) / N * (1 - eps(N))
        close(expected_boot, 1 - 1 / N + np.trace(cov) / (N * N * full_var))
        members = rows(base / f'memberships/{dataset}_bootstrap_multiplicities_v01.csv.gz')
        reps = sorted({r['replicate_id'] for r in members})
        ri, weights = {v: i for i, v in enumerate(reps)}, np.zeros((len(reps), N))
        for r in members:
            weights[ri[r['replicate_id']], ii[r['unit_id']]] = float(r['multiplicity'])
        assert (weights >= 0).all() and (weights.sum(axis=1) == N).all()
        weights /= N
        ratios = np.einsum('bi,ij,bj->b', weights, cov, weights) / full_var
        close(ratios.mean(), boot['variance_ratio_mean'])
        close(np.median(ratios), boot['variance_ratio_median'])
        subset_ref = math.sqrt(n * (1 + (N - 1) * rho) / (N * (1 + (n - 1) * rho)))
        split_ref = rho * math.sqrt(n * m) / math.sqrt((1 + (n - 1) * rho) * (1 + (m - 1) * rho))
        answers.append(dict(dataset=dataset, N=N, T=T, n=n, m=m, rbar=rho, EPS=eps(N),
            subset_model_reference=subset_ref,
            subset_observed_median=float(subset['correlation_median']),
            subset_median_minus_reference=float(subset['correlation_median']) - subset_ref,
            split_model_reference=split_ref,
            split_observed_median=float(split['correlation_median']),
            split_median_minus_reference=float(split['correlation_median']) - split_ref,
            bootstrap_exact_expected_variance_ratio=expected_boot,
            bootstrap_observed_mean=float(boot['variance_ratio_mean']),
            bootstrap_observed_median=float(boot['variance_ratio_median']),
            bootstrap_B=len(reps),
            bootstrap_mcse_of_mean=float(np.std(ratios, ddof=1) / math.sqrt(len(reps))),
            Deff_observed=float(s['effective_dimension']),
            Deff_population_exchangeable=N / (1 + (N - 1) * rho * rho)))
        corr = np.corrcoef(x, rowvar=False)
        pairs = corr[np.triu_indices(N, k=1)]
        vr = float(np.mean((pairs - pairs.mean()) ** 2))
        direct = N * N / float(np.sum(corr * corr))
        identity = N / (1 + (N - 1) * (float(pairs.mean()) ** 2 + vr))
        close(identity, direct)
        close(identity, s['effective_dimension'])
        deff_rows.append(dict(dataset=dataset, N=N, pairs=len(pairs),
            rbar=float(pairs.mean()), vr_divisor_M=vr, Deff_direct=direct, Deff_identity=identity))
        figure_rows.append(dict(dataset=dataset, EPS=float(s['eps']),
            Deff_over_N=float(s['effective_dimension_fraction'])))

    expected_dir = root / 'revision_support/expected_results'
    expected = json.loads((expected_dir / 'eps_baseline_comparison.json').read_text(encoding='utf-8'))
    expected_deff = json.loads((expected_dir / 'deff_dispersion_identity_check.json').read_text(encoding='utf-8'))
    for observed, reference in [(answers, expected['records']), (deff_rows, expected_deff)]:
        assert len(observed) == len(reference)
        for a, b in zip(observed, reference):
            assert a['dataset'] == b['dataset']
            for key, value in a.items():
                if key != 'dataset':
                    close(value, b[key])

    result = dict(status='Verified from the archived complete standardized data and frozen memberships',
        notes=expected['notes'], records=answers, input_sha256=inputs)
    (out / 'eps_baseline_comparison.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    (out / 'deff_dispersion_identity_check.json').write_text(json.dumps(deff_rows, indent=2) + '\n', encoding='utf-8')
    (out / 'figure2_input_values.json').write_text(json.dumps(figure_rows, indent=2) + '\n', encoding='utf-8')

    if not args.skip_figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                             'axes.linewidth': 0.8, 'savefig.facecolor': 'white'})
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        xx = np.arange(3)
        ax.plot(xx, [r['EPS'] for r in figure_rows], linestyle='None', marker='o',
                markersize=6, color='#0072B2', label='Conventional EPS estimate')
        ax.plot(xx, [r['Deff_over_N'] for r in figure_rows], linestyle='None', marker='s',
                markersize=6, color='#D55E00', label=r'$D_{\mathrm{eff}}/N$')
        ax.set_xticks(xx, ['ca533', 'FRAN003', 'Yoshino'])
        ax.set(xlabel='Empirical chronology', ylabel='Coefficient or normalized dimension',
               ylim=(0, 1.02), xlim=(-.25, 2.25))
        ax.legend(frameon=False, loc='lower left', numpoints=1)
        fig.tight_layout()
        fig.savefig(out / 'fig_main03_empirical_structure_reproduced.pdf', bbox_inches='tight')
        plt.close(fig)

    versions = {'python': __import__('platform').python_version(), 'numpy': np.__version__}
    if not args.skip_figure:
        versions['matplotlib'] = matplotlib.__version__
    report = dict(status='passed', datasets=3, bootstrap_draws_per_dataset=2000,
                  new_resampling=False, raw_detrending_rerun=False,
                  comparison='All numerical records match the revision reference within atol=1e-11, rtol=1e-10',
                  figure_generated=not args.skip_figure, versions=versions,
                  figure_note='Visual geometry and numerical values are reproduced; PDF bytes can differ by software version and metadata.')
    (out / 'reproduction_report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
