#!/usr/bin/env python3
"""Validate the finite-N exact INR against quadrature and existing Monte Carlo.

Requirements: Python 3, numpy, scipy. Run with no arguments from any directory.
The existing 48-point, two-million-sample experiment is read, never rerun.
The only new Monte Carlo uses 200,000 full complex vectors per case.
All INR values are linear means; z scores are diagnostics, not proof/tests.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import hyp1f1, roots_legendre


def exact_inr(n: int, r: float, gamma: float, kappa: float = 1.) -> float:
    """E[INR], with gamma=G1/sigma^2 and kappa=Et/Er; N>=2."""
    if n < 2 or not 0 <= r <= 1 or gamma < 0:
        raise ValueError('Require N>=2, 0<=r<=1 and gamma>=0')
    lam = n * gamma
    b = lam * (1 - r)
    return float(kappa * ((n - 1) * lam * hyp1f1(1, n + 1, -lam) / n
                         - (n - 2) * b * hyp1f1(1, n, -b) / (n - 1)))


def positive_integrand(n, a, b, p, q):
    """Unscaled positive Laplace integral on (p,t) in [0,1]^2, q=p*t."""
    m = n - 1
    polynomial = (a*q*(m-b*(p-q))**2 + a*m*q
                  + a*b*((p-q)**2+q*q) + b*p + b*b*p*q)
    return q**(n-2) * np.exp(-a*(1-p)-b*(1-q)) * polynomial


def adaptive_integral(n, r, gamma, cap=100., tolerance=1e-11):
    """Resolve high-SNR concentration using x=1-p, y=p-q.

    The Jacobian is 1/p. Scale U=max(1,N*gamma)*x and
    V=max(1,N*gamma*(1-r))*y, then integrate the triangular domain.
    The cap discards only exponentially suppressed tails; a second cap is
    checked below. This is a positive integral, not a difference of beta terms.
    """
    a, b = n * gamma * r, n * gamma * (1-r)
    sx, sy = max(1., a+b), max(1., b)

    def outer(u):
        x = u / sx
        p = 1-x

        def inner(v):
            y = v / sy
            q = max(0., p-y)
            m = n-1
            polynomial = (a*q*(m-b*y)**2 + a*m*q
                          + a*b*(y*y+q*q) + b*p + b*b*p*q)
            # Evaluate the exponent in scaled coordinates to avoid rounding
            # 1-p to zero when the integral concentrates near p=t=1.
            return (q**(n-2) * np.exp(-(a+b)*x-b*y) * polynomial
                    / (p*sx*sy))

        return quad(inner, 0., min(sy*p, cap), epsabs=tolerance,
                    epsrel=tolerance, limit=200)[0]

    return quad(outer, 0., min(sx, cap), epsabs=tolerance,
                epsrel=tolerance, limit=200)[0]


@lru_cache(maxsize=None)
def gauss_grid(order):
    x, w = roots_legendre(order)
    return (x+1)/2, w/2


def direct_gauss_integral(n, r, gamma):
    """Independent Gauss-Legendre quadrature in the original square domain.

    Increase order with sqrt(N*gamma) to resolve the corner at high SNR.
    Chunk the outer grid to avoid allocating an order-by-order matrix.
    """
    required = max(64., 4*np.sqrt(n*(1+gamma)))
    order = 2**int(np.ceil(np.log2(required)))
    grid, weights = gauss_grid(order)
    a, b = n*gamma*r, n*gamma*(1-r)
    total = 0.
    for start in range(0, order, 64):
        p = grid[start:start+64, None]
        q = p * grid[None, :]
        values = positive_integrand(n, a, b, p, q)
        total += float(weights[start:start+64] @ (values @ weights))
    return total, order


def n2_analytic(gamma):
    """For N=2, E[INR]/kappa=1-(1-exp(-2*gamma))/(2*gamma)."""
    lam = 2*gamma
    if lam == 0:
        return 0.
    if lam < .1:
        # Independent convergent Taylor series avoids cancellation at zero.
        term = lam/2
        value = term
        for j in range(2, 30):
            term *= -lam/(j+1)
            value += term
        return value
    return float(1+np.expm1(-lam)/lam)


def relative_error(value, reference):
    return float(abs(value-reference)/abs(reference))


def deterministic_checks():
    dimensions, overlaps = [2, 4, 16, 64], [0., .4, .999999, 1.]
    n2_errors = [relative_error(exact_inr(2, r, gamma), n2_analytic(gamma))
                 for r in overlaps for gamma in [1e-10, .001, 1., 100., 1e4]]
    assert max(n2_errors) < 1e-12
    for n in dimensions:
        for r in overlaps:
            assert exact_inr(n, r, 0.) == 0.
    low_gamma = 1e-10
    low_errors = [relative_error(exact_inr(n, r, low_gamma)/low_gamma,
                                (1+n*(n-2)*r)/(n-1))
                  for n in dimensions for r in overlaps]
    assert max(low_errors) < 1e-7
    high_gamma = 1e10
    high_errors = [relative_error(exact_inr(n, r, high_gamma),
                                 n-1 if r == 1 else 1.)
                   for n in dimensions for r in [0., .4, 1.]]
    assert max(high_errors) < 1e-7
    return {
        'N2_analytic_cases': len(n2_errors),
        'N2_analytic_max_relative_error': max(n2_errors),
        'zero_channel_cases': len(dimensions)*len(overlaps),
        'low_gamma': low_gamma,
        'low_gamma_slope': '(1+N*(N-2)*r)/(N-1)',
        'low_gamma_max_relative_error': max(low_errors),
        'high_gamma': high_gamma,
        'high_gamma_overlaps': [0., .4, 1.],
        'high_gamma_limit': '1 for fixed r<1; N-1 for r=1 (in units of kappa)',
        'high_gamma_max_relative_error': max(high_errors),
        'note': 'The r<1 high-SNR limit is not uniform as r approaches 1.',
    }


def existing_mc_checks(data, kappa):
    rows = []
    for row in data['rows']:
        exact = exact_inr(row['N'], row['r'], row['gamma'], kappa)
        mean, se = row['exact']['mean'], row['exact']['se']
        rows.append([row['N'], row['r'], row['gamma'], exact, mean, se,
                     (mean-exact)/se])
    z = np.array([row[-1] for row in rows])
    return {
        'source': 'results.json',
        'seed': data['seed'],
        'samples_per_point': data['samples_per_point'],
        'points': len(rows),
        'mean_z': float(np.mean(z)),
        'rms_z': float(np.sqrt(np.mean(z*z))),
        'max_abs_z': float(np.max(abs(z))),
        'worst_case_N_r_gamma': rows[int(np.argmax(abs(z)))][:3],
        'columns': ['N', 'r', 'gamma', 'closed_form', 'mc_mean', 'mc_se', 'z'],
        'rows': rows,
    }


def quadrature_checks():
    rows, cap_errors = [], []
    for n in [2, 4, 16, 64]:
        for r in [0., .4, .999999, 1.]:
            for gamma in [.001, 1., 100., 1e4]:
                exact = exact_inr(n, r, gamma)
                adaptive = adaptive_integral(n, r, gamma)
                gauss, order = direct_gauss_integral(n, r, gamma)
                ea, eg = relative_error(adaptive, exact), relative_error(gauss, exact)
                if gamma == 1e4:
                    tight = adaptive_integral(n, r, gamma, cap=120., tolerance=1e-12)
                    cap_errors.append(relative_error(tight, adaptive))
                assert ea < 2e-8, (n, r, gamma, 'adaptive', ea)
                assert eg < 2e-6, (n, r, gamma, 'direct Gauss', eg)
                rows.append([n, r, gamma, exact, ea, eg, order])
        print(f'Quadrature checked N={n}', flush=True)
    assert max(cap_errors) < 2e-8
    return {
        'points': len(rows),
        'adaptive_cap': 100.,
        'adaptive_tolerance': 1e-11,
        'max_adaptive_relative_error': max(row[4] for row in rows),
        'max_direct_gauss_relative_error': max(row[5] for row in rows),
        'high_gamma_cap120_tolerance1e12_max_relative_change': max(cap_errors),
        'columns': ['N', 'r', 'gamma', 'closed_form_over_kappa',
                    'adaptive_relative_error', 'direct_gauss_relative_error',
                    'direct_gauss_order'],
        'rows': rows,
    }


def full_vector_mc(n, r, gamma, kappa, rng, samples=200_000):
    """Direct normalized projection, dense complex channels, no reduction."""
    unitary, _ = np.linalg.qr(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
    h0 = np.zeros(n, dtype=complex)
    h0[0] = 1.
    h1 = np.zeros(n, dtype=complex)
    h1[:2] = np.sqrt(n*gamma)*np.array([np.sqrt(r)*np.exp(.73j), np.sqrt(1-r)])
    rotated_h0, rotated_h1 = unitary @ h0, unitary @ h1
    total, total_sq, rotation_error = 0., 0., 0.

    def values(y, h0, h1):
        q = np.sum(abs(y)**2, axis=1)
        projected = h0-y*((y.conj() @ h0)/q)[:, None]
        beam = projected / np.linalg.norm(projected, axis=1)[:, None]
        return kappa*abs(beam.conj() @ h1)**2

    for start in range(0, samples, 20_000):
        size = min(20_000, samples-start)
        noise = (rng.normal(size=(size, n))+1j*rng.normal(size=(size, n)))/np.sqrt(2)
        y = h1+noise
        original = values(y, h0, h1)
        rotated = values(y @ unitary.T, rotated_h0, rotated_h1)
        rotation_error = max(rotation_error, float(np.max(abs(original-rotated))))
        total += float(np.sum(rotated))
        total_sq += float(np.sum(rotated**2))
    mean = total/samples
    se = np.sqrt((total_sq-samples*mean**2)/(samples-1)/samples)
    exact = exact_inr(n, r, gamma, kappa)
    assert rotation_error < 1e-11
    return {'N': n, 'r': r, 'gamma': gamma, 'samples': samples,
            'closed_form': exact, 'mc_mean': mean, 'mc_se': float(se),
            'z': float((mean-exact)/se),
            'unitary_invariance_max_absolute_error': rotation_error}


def main():
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=directory/'results.json')
    parser.add_argument('--output', type=Path, default=directory/'exact_results.json')
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    kappa = data['Et']/data['Er']
    result = {
        'formula': 'kappa*((N-1)*L*beta_N(L)-(N-2)*B*beta_(N-1)(B))',
        'definitions': 'L=N*gamma; B=L*(1-r); beta_n(x)=hyp1f1(1,n+1,-x)/n',
        'kappa': kappa,
        'deterministic': deterministic_checks(),
        'existing_mc': existing_mc_checks(data, kappa),
        'quadrature': quadrature_checks(),
    }
    seed = 20260923
    rng = np.random.default_rng(seed)
    result['fresh_full_vector_mc'] = {
        'seed': seed,
        'rows': [full_vector_mc(n, .4, 1., kappa, rng) for n in [2, 16]],
        'note': 'Monte Carlo z scores are diagnostics, not fixed pass/fail thresholds.',
    }
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({
        'existing_mc_points': result['existing_mc']['points'],
        'existing_mc_max_abs_z': result['existing_mc']['max_abs_z'],
        'max_adaptive_relative_error': result['quadrature']['max_adaptive_relative_error'],
        'max_direct_gauss_relative_error': result['quadrature']['max_direct_gauss_relative_error'],
        'fresh_mc_z': [row['z'] for row in result['fresh_full_vector_mc']['rows']],
        'output': args.output.name,
    }, indent=2))


if __name__ == '__main__':
    main()
