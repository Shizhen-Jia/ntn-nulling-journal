#!/usr/bin/env python3
"""Reproduce the INR equation 27/28 denominator-ablation experiment.

Requirements: Python 3, numpy, scipy, matplotlib.

Examples:
    python verify_inr_approximations.py --plot-only
    python verify_inr_approximations.py --samples 2000000 --seed 20260922
    python verify_inr_approximations.py --samples 2000 --output /tmp/inr_smoke.json

The default run uses 2,000,000 independent samples per point, with exactly the
ULA geometry in verify_equations.ipynb cell 1.  All reported INR means are
arithmetic means in linear units, converted to dB only for plotting.  The
noise variance N0 is fixed, not estimated anew in each realization.

A unitary basis change reduces the full N-dimensional Gaussian experiment to
two complex Gaussian coordinates and a Gamma(N-2, sigma^2) residual norm.  This
is an exact distributional reduction, not a large-array approximation:
    h0 = e0, h1 = sqrt(N*G1) * (sqrt(r)*e0 + sqrt(1-r)*e1),
    v ~ CN(0, sigma^2 I), sigma^2 = N0/Er, hhat = h1 + v.
Scaling h0 has no effect on its normalized nulling beamformer.  Gaussian
rotational invariance means the resulting law depends on the original fixed
channels only through N, G1 and r = |rho|^2.

"Frozen q" deliberately replaces the internal q = ||hhat||^2 by E[q] in the
leakage expression.  It is a diagnostic surrogate, NOT a feasible normalized
nulling beamformer.  Eq. (28) is labeled an approximation, not an exact answer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import hyp1f1

ER, ET, N0 = 100.0, 1.0, 0.1
SIGMA2 = N0 / ER
DEFAULT_SAMPLES = 2_000_000
DEFAULT_SEED = 20260922


def statistics(values: np.ndarray) -> dict[str, float]:
    """Standard error of the linear sample mean, not scatter in dB."""
    return {
        'mean': float(values.mean()),
        'se': float(values.std(ddof=1) / np.sqrt(len(values))),
        'db': float(10 * np.log10(values.mean())),
    }


def exact_expected_rho2(n: int, r: float, gamma: float) -> float:
    """Exact E|rhohat|^2 for fixed channels and circular Gaussian noise.

    Set lambda = N*gamma and K ~ Poisson(lambda).  In coordinates parallel
    and perpendicular to the true h1 direction, the mean normalized channel
    projector has eigenvalues 1-(N-1)b and b, respectively, where
        b = E[1/(N+K)] = integral_0^1 t^(N-1) exp(lambda*(t-1)) dt
          = hyp1f1(1, N+1, -lambda)/N.
    The hypergeometric form avoids an overflowing exp(lambda) factor and
    needs no Monte Carlo estimate of the denominator being frozen.
    """
    b = hyp1f1(1, n + 1, -n * gamma) / n
    return float(r * (1 - (n - 1) * b) + (1 - r) * b)


def theory(n: int, r: float, gamma, sign: float):
    return ((ET / ER) * gamma * (n * r + gamma * (1 + sign * r)) /
            ((1 + gamma) * (1 + gamma * (1 - r))))


def run_case(n: int, r: float, gamma: float, samples: int,
             rng: np.random.Generator) -> dict:
    channel_power = n * SIGMA2 * gamma
    c = np.sqrt(channel_power * r)
    d = np.sqrt(channel_power * (1 - r))
    v0 = np.sqrt(SIGMA2 / 2) * (
        rng.standard_normal(samples) + 1j * rng.standard_normal(samples))
    v1 = np.sqrt(SIGMA2 / 2) * (
        rng.standard_normal(samples) + 1j * rng.standard_normal(samples))
    rest = (rng.gamma(shape=n - 2, scale=SIGMA2, size=samples)
            if n > 2 else np.zeros(samples))
    z0, z1 = c + v0, d + v1
    q = np.abs(z0)**2 + np.abs(z1)**2 + rest
    # This positive-sum form avoids cancellation in 1-|z0|^2/q near r=1.
    a = (np.abs(z1)**2 + rest) / q
    hhat_h = np.conj(z0) * c + np.conj(z1) * d
    expected_rho2 = exact_expected_rho2(n, r, gamma)
    expected_a = 1 - expected_rho2
    expected_q = channel_power + n * SIGMA2
    raw = c - z0 * hhat_h / q
    if gamma >= 100 and r < 1:
        # Exact nulling gives w^H h1 = -w^H v.  This identity avoids
        # subtracting large channel terms to obtain tiny high-SNR leakage.
        hhat_v = np.conj(z0) * v0 + np.conj(z1) * v1 + rest
        raw = -(v0 - z0 * hhat_v / q)
    elif r == 1:
        raw = c * a
    raw_frozen = c - z0 * hhat_h / expected_q
    factor = ET / N0
    row = {
        'N': n, 'r': float(r), 'gamma': float(gamma),
        'G1': SIGMA2 * gamma, 'exact_erho': expected_rho2,
    }
    for name, values in [
        ('exact', factor * np.abs(raw)**2 / a),
        ('outer_frozen', factor * np.abs(raw)**2 / expected_a),
        ('inner_frozen', factor * np.abs(raw_frozen)**2 / a),
        ('both_frozen', factor * np.abs(raw_frozen)**2 / expected_a),
    ]:
        row[name] = statistics(values)
    row['eq27'] = float(theory(n, r, gamma, 1))
    row['eq28'] = float(theory(n, r, gamma, -1))
    # Exact second moment for this deliberately doubly-frozen surrogate,
    # retaining the quadratic-noise variance and finite-N mean correction.
    row['both_frozen_exact_mean'] = float(factor / expected_a * (
        c*c * ((n - 1)*SIGMA2 / expected_q)**2 +
        SIGMA2 * channel_power**2 * (1 + r) / expected_q**2 +
        SIGMA2**2 * channel_power / expected_q**2))
    row['relative_eq28_error'] = row['eq28'] / row['exact']['mean'] - 1
    print(f'N={n:2d}, r={r:.8f}, gamma={gamma:g}: '
          f"exact={row['exact']['mean']:.8g} +/- {row['exact']['se']:.3g}, "
          f"eq27={row['eq27']:.8g}, eq28={row['eq28']:.8g}", flush=True)
    return row


def run_experiment(samples: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    rows = []
    for n in [16, 32, 48, 64]:
        h0 = np.exp(1j * np.pi * np.arange(n))
        h1 = np.exp(1j * np.pi * np.arange(n) * np.cos(np.deg2rad(15)))
        r = float(abs(np.vdot(h1, h0))**2 / n**2)
        for gamma in [.001, .1, 1, 10, 100, 1e4]:
            rows.append(run_case(n, r, gamma, samples, rng))
    for n in [4, 16, 64]:
        for r in [0., 1.]:
            for gamma in [.001, 1, 100, 1e4]:
                rows.append(run_case(n, r, gamma, samples, rng))
    return {
        'seed': seed, 'samples_per_point': samples,
        'Er': ER, 'Et': ET, 'N0': N0, 'rows': rows,
    }


def check_identities(data: dict) -> dict:
    """Deterministic algebra check and MC-vs-analytic frozen-model check."""
    rng = np.random.default_rng(901)
    n, r, gamma, samples = 16, .778974957047383, 10., 10000
    power = n * SIGMA2 * gamma
    zeta = np.sqrt(SIGMA2 / 2) * (
        rng.normal(size=samples) + 1j * rng.normal(size=samples))
    x = np.sqrt(SIGMA2 / 2) * (
        rng.normal(size=samples) + 1j * rng.normal(size=samples))
    q_perp = np.abs(x)**2 + rng.gamma(n - 2, SIGMA2, size=samples)
    y_parallel = np.sqrt(power) + zeta
    q = np.abs(y_parallel)**2 + q_perp
    h0_hat = np.sqrt(r)*y_parallel + np.sqrt(1-r)*x
    direct = np.sqrt(power*r) - h0_hat*np.conj(y_parallel)*np.sqrt(power)/q
    reduced = np.sqrt(power) * (
        np.sqrt(r)*q_perp - np.sqrt(1-r)*x*np.conj(y_parallel)) / q
    identity_error = float(np.max(np.abs(direct - reduced)))
    if identity_error > 1e-12:
        raise AssertionError(f'Exact numerator identity failed: {identity_error}')
    z_scores = [abs(row['both_frozen']['mean'] - row['both_frozen_exact_mean']) /
                row['both_frozen']['se'] for row in data['rows']]
    return {
        'numerator_identity_max_absolute_error': identity_error,
        'both_frozen_max_abs_mean_error_in_SE': float(max(z_scores)),
        'points_checked': len(z_scores),
        'note': 'The MC z score is a diagnostic, not a proof or a fixed test threshold.',
    }


def make_plot(data: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6), sharey=True)
    gamma = np.logspace(-5, 6, 1000)
    for ax, n in zip(axes.flat, [16, 32, 48, 64]):
        rows = [x for x in data['rows'] if x['N'] == n and 0 < x['r'] < 1]
        rows.sort(key=lambda row: row['gamma'])
        r = rows[0]['r']
        ax.semilogx(SIGMA2*gamma, 10*np.log10(theory(n, r, gamma, 1)),
                    '--', color='#2166ac', lw=1.8, label='Approx. eq. (27), plus')
        ax.semilogx(SIGMA2*gamma, 10*np.log10(theory(n, r, gamma, -1)),
                    '-', color='#d73027', lw=1.8, label='Approx. eq. (28), minus')
        for key, marker, color, label, zorder in [
            ('inner_frozen', '^', '#762a83', 'Frozen q: diagnostic surrogate', 4),
            ('outer_frozen', 's', '#1b9e77', 'Frozen outer normalization only', 5),
            ('exact', 'o', '#111111', 'Exact beamformer MC (95% CI)', 6),
        ]:
            xs = np.array([row['G1'] for row in rows])
            mean = np.array([row[key]['mean'] for row in rows])
            se = np.array([row[key]['se'] for row in rows])
            ys = 10*np.log10(mean)
            lower = np.maximum(mean-1.96*se, mean*1e-6)
            errors = np.vstack([ys-10*np.log10(lower),
                                10*np.log10(mean+1.96*se)-ys])
            ax.errorbar(xs, ys, yerr=errors, linestyle='none', marker=marker,
                        markersize=6, markerfacecolor='none', markeredgewidth=1.2,
                        capsize=3, color=color, label=label, zorder=zorder)
        ax.set_title(f'N = {n},  |rho|^2 = {r:.6f}')
        ax.set_xlabel(r'$G_1$  ($\gamma_r=1000G_1$)')
        ax.set_ylabel('Mean INR [dB]')
        ax.set_xlim(SIGMA2*gamma[0], SIGMA2*gamma[-1])
        ax.set_ylim(-80, -7)
        ax.grid(True, alpha=.35)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .97),
               ncol=2, frameon=False, fontsize=10)
    fig.suptitle('INR: denominator ablation and equations (27)/(28)', y=.995, fontsize=15)
    count = data['samples_per_point']
    fig.text(.5, .012,
             f'{count:,} independent samples per point; seed {data["seed"]}. '
             'Confidence intervals may be smaller than markers.\n'
             'Frozen q is a leakage diagnostic, not a feasible normalized nulling beamformer.',
             ha='center', va='bottom', fontsize=9)
    fig.tight_layout(rect=(0, .065, 1, .875))
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f'Wrote plot: {path}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--samples', type=int, default=DEFAULT_SAMPLES)
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().with_name('results.json'))
    parser.add_argument('--plot-only', action='store_true',
                        help='Read --output and regenerate plot; do not rerun Monte Carlo.')
    args = parser.parse_args()
    if args.samples < 2:
        parser.error('--samples must be at least 2')
    output = args.output.expanduser().resolve()
    if args.plot_only:
        with output.open() as f:
            data = json.load(f)
    else:
        data = run_experiment(args.samples, args.seed)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('w') as f:
            json.dump(data, f, indent=2)
            f.write('\n')
        print(f'Wrote results: {output}')
    print('Validation summary:', json.dumps(check_identities(data), indent=2))
    # Default output uses comparison.png. A custom output stem keeps smoke
    # runs or alternative experiments from overwriting that deliverable.
    plot_path = output.with_name('comparison.png' if output.name == 'results.json'
                                 else output.stem + '_comparison.png')
    make_plot(data, plot_path)


if __name__ == '__main__':
    main()
