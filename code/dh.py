#!/usr/bin/env python3
"""
dh.py — Davenport–Heilbronn gates G1/G2/G3 (prereg v2.0, §4) with LOGGED
results and the 45-minute wall-clock box, plus t=0 zero location
(window seeding reuses track.seed_zeros; targeted off-axis search below).

Gate thresholds (prereg §4, literal):
  G1: |Λ_DH(s) - Λ_DH(1-s)| < 1e-20 at 5 random s.
  G2: |Im Λ_DH(1/2+it)| < 1e-20 at 5 random t.
  G3: |∫Φ_DH cos - Λ_DH((1+iz)/2)| < 1e-25 up to ONE fitted overall scale,
      at 5 random z.
On timeout or any gate failure: DH is DROPPED (logged red flag); the study
proceeds as ζ vs (b), (c), (e).

NOTE on prereg §4 vs the frozen ht_flow.py: the preregistration writes the DH
Φ sum as n=1..20; ht_flow.phi_dh_wide (frozen, do-not-modify) sums n=1..30.
More terms can only help the G3 cross-check, and G3 itself is the arbiter, so
the 30-term sum is used and the discrepancy is logged, not silently "fixed".
"""
import json
import os
import random
import time
from mpmath import mp

mp.dps = 50

from ht_flow import HeatFlow, phi_dh_wide, dh_kappa, dh_f, dh_Lambda

DH_BOX_SECONDS = 45 * 60  # prereg §4 / §10

# prereg §4: J=256 on [0,3], dps=50
DH_J, DH_UMAX, DH_DPS = 256, 3.0, 50


class BoxExpired(Exception):
    pass


def _check_box(t0):
    if time.time() - t0 > DH_BOX_SECONDS:
        raise BoxExpired("45-minute DH box expired")


def run_gates(log_path=None, seed=7):
    """Run G1/G2/G3. Returns (decision, log) where decision is 'IN' or
    'DROPPED'. The log dict is written to log_path as JSON (prereg: gates
    must be LOGGED). All randomness is seeded (seed=7) for reproducibility."""
    t0 = time.time()
    rng = random.Random(seed)
    log = {'gates': {}, 'seed': seed, 'box_seconds': DH_BOX_SECONDS,
           'note': 'prereg §4: DH Φ sum written n=1..20; frozen ht_flow.py '
                   'implements n=1..30 (more terms; G3 is the arbiter)'}
    try:
        kappa = dh_kappa()
        log['kappa'] = str(kappa)
        log['kappa_expect'] = '0.284079043840...'
        _check_box(t0)

        # G1: functional equation
        worst = mp.mpf(0)
        pts = []
        for _ in range(5):
            s = mp.mpc(rng.uniform(-1, 2), rng.uniform(1, 50))
            r = abs(dh_Lambda(s, kappa) - dh_Lambda(1 - s, kappa))
            worst = max(worst, r)
            pts.append([float(mp.re(s)), float(mp.im(s)), float(r)])
            _check_box(t0)
        g1 = {'worst': float(worst), 'threshold': 1e-20,
              'pass': bool(worst < 1e-20), 'points': pts}
        log['gates']['G1'] = g1

        # G2: reality on the critical line
        worst = mp.mpf(0)
        pts = []
        for _ in range(5):
            tt = mp.mpf(rng.uniform(1, 100))
            r = abs(mp.im(dh_Lambda(mp.mpf('0.5') + mp.j * tt, kappa)))
            worst = max(worst, r)
            pts.append([float(tt), float(r)])
            _check_box(t0)
        g2 = {'worst': float(worst), 'threshold': 1e-20,
              'pass': bool(worst < 1e-20), 'points': pts}
        log['gates']['G2'] = g2

        # G3: quadrature cross-check up to one fitted overall scale
        hf = HeatFlow(lambda u: phi_dh_wide(u, kappa),
                      J=DH_J, umax=DH_UMAX, dps=DH_DPS)
        log['dh_quadrature_J_kept'] = hf.J
        z_tests = [mp.mpf(zr) for zr in
                   [str(rng.uniform(3, 160)) for _ in range(5)]]
        ratios = []
        for z in z_tests:
            q = hf.H(z, 0)
            lam = dh_Lambda((1 + mp.j * z) / 2, kappa)
            ratios.append(q / lam)
            _check_box(t0)
        r0 = ratios[0]
        spread = max(abs(r / r0 - 1) for r in ratios)
        g3 = {'spread': float(spread), 'threshold': 1e-25,
              'pass': bool(spread < 1e-25),
              'z_tests': [float(z) for z in z_tests],
              'scale_ratios': [str(r) for r in ratios]}
        log['gates']['G3'] = g3

        # PDE sanity (not a gate; logged for the record)
        try:
            log['pde_residual_37_0p1'] = float(
                hf.pde_residual(mp.mpf('37'), mp.mpf('0.1')))
        except Exception as e:  # noqa: BLE001
            log['pde_residual_37_0p1'] = 'error: %s' % e
        _check_box(t0)

        all_pass = g1['pass'] and g2['pass'] and g3['pass']
        log['decision'] = 'IN' if all_pass else 'DROPPED'
        log['elapsed_s'] = time.time() - t0
    except BoxExpired as e:
        log['decision'] = 'DROPPED'
        log['box_expired'] = str(e)
        log['elapsed_s'] = time.time() - t0
    if log_path:
        with open(log_path, 'w') as f:
            json.dump(log, f, indent=1)
    return log['decision'], log


def make_dh_flow(kappa):
    """HeatFlow for the DH Φ (prereg §4 parameters). t=0 use only until
    gates pass."""
    return HeatFlow(lambda u: phi_dh_wide(u, kappa),
                    J=DH_J, umax=DH_UMAX, dps=DH_DPS)


def find_offaxis_zeros(kappa, log_path=None, box_seconds=DH_BOX_SECONDS):
    """Targeted Newton search for genuinely off-critical-line zeros of f(s).

    Scans a grid of (sigma, t) seeds; keeps converged zeros with
    |Re(s) - 1/2| > 1e-3 AND |f(s)| < 1e-20. v1's "spira" seed polished to
    Im ~ 1e-81 (on the real axis) — this scan requires genuine off-axis.
    Returns (found_list, log); found entries are wide-frame z = 2*Im(s) -
    i*(2*Re(s)-1). If empty, the off-axis pair is NOT located (prereg R3:
    logged, DH comparison rests on the window signature).
    """
    t0 = time.time()
    from ht_flow import dh_f as _dh_f
    found = []
    log = {'seeds_tried': 0, 'found': []}
    for sigma in [0.3, 0.4, 0.6, 0.7, 0.8]:
        for t_guess in [20.0, 40.0, 60.0, 85.0, 110.0, 140.0]:
            if time.time() - t0 > box_seconds:
                log['box_expired'] = True
                break
            log['seeds_tried'] += 1
            s = mp.mpc(sigma, t_guess)
            try:
                for _ in range(60):
                    f = _dh_f(s, kappa)
                    h = mp.mpf('1e-8')
                    fp = (_dh_f(s + h, kappa) - _dh_f(s - h, kappa)) / (2 * h)
                    if fp == 0:
                        break
                    ds = f / fp
                    s -= ds
                    if abs(ds) < mp.mpf('1e-30'):
                        break
                if abs(_dh_f(s, kappa)) < mp.mpf('1e-20') and \
                        abs(mp.re(s) - mp.mpf('0.5')) > mp.mpf('1e-3'):
                    z = 2 * mp.im(s) - mp.j * (2 * mp.re(s) - 1)
                    if not any(abs(z - w) < mp.mpf('1e-3') for w in found):
                        found.append(z)
                        log['found'].append(
                            {'s': [float(mp.re(s)), float(mp.im(s))],
                             'z': [float(mp.re(z)), float(mp.im(z))]})
            except Exception:  # noqa: BLE001 - one bad seed must not kill the scan
                continue
    log['n_found'] = len(found)
    if log_path:
        with open(log_path, 'w') as f:
            json.dump(log, f, indent=1)
    return found, log


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    data = os.path.join(here, '..', 'data')
    os.makedirs(data, exist_ok=True)
    decision, log = run_gates(os.path.join(data, 'dh_gates.json'))
    print("DH gates decision:", decision)
    for g, r in log['gates'].items():
        print("  %s: %s" % (g, 'PASS' if r['pass'] else 'FAIL'))
