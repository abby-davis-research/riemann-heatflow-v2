#!/usr/bin/env python3
"""run_gates.py — prereg v2.1 §7 gates V1, V2, V2b, V3, V4, V6 + §8(i).

V2b (new in v2.1): seeding completeness — the v2.0 failure mode.
V5 (seed independence) needs the Phase-1 down-leg and runs at the end of
run_phase1.py; it is recorded as PENDING here. Phase 2 is barred until the
gate log shows every gate passing. DH case: DROPPED per prereg v2.1 §4
(carried from v2.0; G3 spread 3.2e-07 vs 1e-25, separate DH-quadrature
issue, out of scope for the v2.1 instrument repair).

Exit code 0: all runnable gates pass. Non-zero: a gate failed -> STOP.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mpmath import mp  # noqa: E402
import run_common as C  # noqa: E402
import analyze  # noqa: E402
import track  # noqa: E402
from synthetic import SyntheticFlow  # noqa: E402
# (dh not imported in v2.1: DH case DROPPED per prereg §4, not re-gated)


def main():
    t_start = time.time()
    print('=== v2.1 gates: V1, V2, V2b, V3, V4, V6, 8(i) ===', flush=True)
    f256, f512 = C.build_zeta_flows()

    # t=0 seeding (prereg §5)
    seeds = track.label_zeros(track.seed_zeros(f256, 300.0))
    C.save_seeds(os.path.join(C.DATADIR, 'seeds_zeta.json'), seeds)
    print('zeta seeds: %d zeros in |z|<=300' % len(seeds), flush=True)
    # reference: 52 positive zeta zeros with gamma <= 150 -> 104 in |z|<=300
    print('  (reference: 104 expected from mpmath zetazero; V2b gates this)',
          flush=True)

    results = {}

    print('-- V1 (quadrature J=256 -> J=512) --', flush=True)
    results['V1'] = analyze.gate_V1(f256, f512)
    print('   worst dH = %.2e (< 1e-60): %s'
          % (results['V1']['worst'], results['V1']['pass']), flush=True)

    print('-- V2 (t=0 zeros, first five) --', flush=True)
    results['V2'] = analyze.gate_V2(f256)
    print('   worst |H0| = %.2e (< 1e-30), worst rel = %.2e (< 1e-14): %s'
          % (results['V2']['worst_H0'], results['V2']['worst_rel'],
             results['V2']['pass']), flush=True)

    print('-- V2b (seeding completeness, NEW in v2.1) --', flush=True)
    results['V2b'] = analyze.gate_V2b(f256)
    print('   seeded %d / %d refs, worst disp = %.2e, worst miss = %.2e: %s'
          % (results['V2b']['n_seeded'], results['V2b']['n_reference'],
             results['V2b']['worst_displacement'],
             results['V2b']['worst_reference_miss'], results['V2b']['pass']),
          flush=True)
    if not results['V2b']['pass']:
        C.add_red_flag('V2b-FAIL',
                       'seeding incomplete: %d/%d, disp %.2e, miss %.2e'
                       % (results['V2b']['n_seeded'],
                          results['V2b']['n_reference'],
                          results['V2b']['worst_displacement'],
                          results['V2b']['worst_reference_miss']))

    print('-- V3 (PDE residual, analytic dH/dt) --', flush=True)
    results['V3'] = analyze.gate_V3(f256, 'zeta')
    print('   worst = %.2e (< 1e-60): %s'
          % (results['V3']['worst'], results['V3']['pass']), flush=True)

    print('-- V4 (synthetic M=48 vs M=96) --', flush=True)
    z1 = 2 * mp.zetazero(1).imag
    syn48 = SyntheticFlow(f256, z1, 1.0, M=48)
    results['V4'] = analyze.gate_V4(syn48, syn48.w)
    print('   worst dH = %.2e, construction = %.2e, pde = %.2e: %s'
          % (results['V4']['worst_Mdiff'], results['V4']['construction'],
             results['V4']['pde'], results['V4']['pass']), flush=True)

    # First t!=0 evaluation happens inside V6 -> the 48h study box starts now.
    C.stamp_study_start()
    print('-- V6 (detector null on zeta (0.2, 0.5]) --', flush=True)
    results['V6'] = analyze.gate_V6(f256, seeds)
    print('   triggers = %d (require 0): %s'
          % (results['V6']['n_triggers'], results['V6']['pass']), flush=True)
    if not results['V6']['pass']:
        C.add_red_flag('V6-FAIL',
                       'detector fired %d time(s) in the rigorous null window'
                       % results['V6']['n_triggers'])

    # §8(i): N_off(0) == 0 for zeta
    n_off0 = sum(1 for z in seeds.values()
                 if abs(mp.im(mp.mpc(z))) > 1e-6)
    results['CAL_i'] = {'gate': 'CAL-(i)', 'N_off_0': n_off0,
                        'pass': n_off0 == 0}
    print('-- CAL-(i): N_off(0) = %d (require 0): %s'
          % (n_off0, n_off0 == 0), flush=True)

    # V5: pending until Phase 1 produces the down-leg
    results['V5'] = {'gate': 'V5', 'pass': None, 'status': 'pending',
                     'note': 'runs at end of run_phase1.py against the '
                             'zeta down-leg'}

    # DH case: DROPPED per prereg v2.1 §4 (carried from v2.0; G3 spread
    # 3.2e-07 vs 1e-25 is a separate DH-quadrature issue, out of scope for
    # the v2.1 instrument repair). Comparison proceeds as zeta vs
    # synth-d1, synth-d01, null-window. Not re-run here.
    results['DH'] = {'gate': 'DH-case', 'decision': 'DROPPED',
                     'pass': True,
                     'note': 'DROPPED per prereg v2.1 §4 (carried from v2.0); '
                             'comparison proceeds as zeta vs synth-d1, '
                             'synth-d01, null-window'}
    print('   DH decision: DROPPED (per prereg v2.1 §4)', flush=True)
    C.add_red_flag('DH-DROPPED',
                   'DH case dropped per prereg v2.1 §4 (carried from v2.0)')

    ok = analyze.log_gates(os.path.join(C.DATADIR, 'gate_log.json'), results)
    print('gate log -> data/gate_log.json; V1,V2,V2b,V3,V4,V6,CAL-(i) pass: %s'
          % all(results[g]['pass'] for g in ('V1', 'V2', 'V2b', 'V3', 'V4',
                                            'V6', 'CAL_i')), flush=True)
    print('elapsed: %.1f min' % ((time.time() - t_start) / 60), flush=True)
    if not all(results[g]['pass']
               for g in ('V1', 'V2', 'V2b', 'V3', 'V4', 'V6', 'CAL_i')):
        print('GATE FAILURE -> STOP. Phase 2 barred.', flush=True)
        return 1
    print('Gates pass. V5 pending Phase 1. Proceed to run_phase1.py.', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
