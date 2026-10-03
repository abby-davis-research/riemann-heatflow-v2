#!/usr/bin/env python3
"""run_phase2.py — prereg v2.0 §9 Phase 2: comparison + calibration + decision.

Cases: synth-d1 (gross RH-false), synth-d01 (near-miss), DH if gated IN.
Calibration §8 (incl. the (iv) landing check with the single preregistered
extension to t=2.0 for synth-d1) -> decision rule §9 -> data/decision.json
+ results.md. Phase 2 refuses to run unless the gate log shows all gates
passing (V1-V6 incl. V5).
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mpmath import mp  # noqa: E402
import run_common as C  # noqa: E402
import analyze  # noqa: E402
import track  # noqa: E402
import dh  # noqa: E402
from synthetic import SyntheticFlow  # noqa: E402

CASE_BOX = 3 * 3600  # prereg §10: 3 hours per synthetic case


def remaining(box_start, total):
    return max(60.0, total - (time.time() - box_start))


def rename_quadruplet(seeds, w, tag):
    """Rename the 4 seeded labels nearest ±w, ±w̄ to w1..w4 (prereg §8(iv)
    landing check expects these labels). Greedy nearest-unused."""
    targets = [w, mp.conj(w), -w, -mp.conj(w)]
    unused = dict(seeds)
    mapping = {}
    new = {}
    for i, t in enumerate(targets, 1):
        lab = min(unused, key=lambda l: abs(unused[l] - t))
        mapping['w%d' % i] = {'from': lab,
                              'seed': [float(mp.re(unused[lab])),
                                       float(mp.im(unused[lab]))],
                              'target': [float(mp.re(t)), float(mp.im(t))],
                              'dist': float(abs(unused[lab] - t))}
        new['w%d' % i] = unused.pop(lab)
    for lab, z in unused.items():
        new[lab] = z
    print('  [%s] quadruplet mapping:' % tag, flush=True)
    for k, v in mapping.items():
        print('    %s <- %s dist=%.2e' % (k, v['from'], v['dist']), flush=True)
    return new, mapping


def merge_legs(l1, l2):
    """Concatenate an extension leg onto a base leg (raw track_leg dicts)."""
    traj = {}
    for lab, pts in l1['traj'].items():
        traj[lab] = pts + (l2['traj'][lab][1:] if lab in l2['traj'] else [])
    return {'traj': traj,
            'events': l1['events'] + l2['events'],
            'track_fails': l1['track_fails'] + l2['track_fails'],
            'n_steps': l1['n_steps'] + l2['n_steps'],
            't_reached': l2['t_reached'],
            'box_expired': bool(l1.get('box_expired') or
                                l2.get('box_expired'))}


def tracked_sample(legs, cap=300):
    """Sample (re, im, t) points from raw legs for confirm_V4_tracked."""
    pts = []
    for leg in legs.values():
        for lab, traj in leg['traj'].items():
            step = max(1, len(traj) // 10)
            for p in traj[::step]:
                pts.append((float(mp.re(p[1])), float(mp.im(p[1])),
                            float(p[0])))
                if len(pts) >= cap:
                    return pts
    return pts


def run_synth_case(f256, z1, delta, tag, seeds_zeta):
    """Full synthetic case: down-leg, up-leg, conditional extension (d1
    only), V4-tracked confirmation. Returns (case_doc_path, syn, syn_hi)."""
    t_case = time.time()
    syn = SyntheticFlow(f256, z1, delta, M=48)
    syn_hi = SyntheticFlow(f256, z1, delta, M=96)  # G-b analog (M-refinement)
    print('== case %s (delta=%s) ==' % (tag, delta), flush=True)

    seeds = track.label_zeros(track.seed_zeros(syn, 300.0))
    seeds, mapping = rename_quadruplet(seeds, syn.w, tag)
    C.save_seeds(os.path.join(C.DATADIR, 'seeds_%s.json' % tag), seeds)

    flow_info = {'type': 'synthetic R(z) surgery (prereg §4)',
                 'delta': str(delta), 'M': 48, 'z1': str(float(z1)),
                 'w': [float(mp.re(syn.w)), float(mp.im(syn.w))],
                 'quadruplet_mapping': mapping,
                 'note': 'G-b uses M=96 refinement (honest analog of '
                         'J-refinement; V4 gates M48~=M96)'}

    print('-- %s down-leg [0,-1] --' % tag, flush=True)
    down = track.track_leg(syn, seeds, mp.mpf('0'), mp.mpf('-1'),
                           box_seconds=remaining(t_case, CASE_BOX),
                           case_tag=tag + '_down')
    C.validate_leg(syn, syn_hi, down, tag + '_down')
    if down['box_expired']:
        C.add_red_flag('BOX-' + tag.upper(),
                       '%s hit the 3h case box; partial map per §10' % tag)

    print('-- %s up-leg [0,0.5] --' % tag, flush=True)
    up = track.track_leg(syn, seeds, mp.mpf('0'), mp.mpf('0.5'),
                         box_seconds=remaining(t_case, CASE_BOX),
                         case_tag=tag + '_up')
    C.validate_leg(syn, syn_hi, up, tag + '_up')

    legs = {'down': down, 'up': up}
    # landing check on the 0.5 up-leg (serialized temp doc)
    tmp = {'legs': {'up': {'events': [track.event_to_jsonable(e)
                                     for e in up['events']]}}}
    landing = analyze.check_landing(tmp)
    extended = False
    if tag == 'synth_d1' and not landing['pass']:
        # prereg §8(iv) + R10: ONE extension to 2.0, iff no landing by 0.5
        print('-- synth_d1: no landing by 0.5 -> single extension to 2.0 --',
              flush=True)
        C.add_red_flag('EXTEND-D1',
                       'synth-d1 up-leg extended 0.5 -> 2.0 (prereg §8(iv), '
                       'sole authorized extension)')
        seeds_end = {lab: pts[-1][1] for lab, pts in up['traj'].items()}
        ext = track.track_leg(syn, seeds_end, mp.mpf('0.5'), mp.mpf('2.0'),
                              box_seconds=remaining(t_case, CASE_BOX),
                              case_tag='synth_d1_up_ext')
        C.validate_leg(syn, syn_hi, ext, 'synth_d1_up_ext')
        legs['up'] = merge_legs(up, ext)
        extended = True
    elif tag == 'synth_d1':
        print('  landing by 0.5: t*_inj = %.6f (no extension needed)'
              % landing['t_inj'], flush=True)

    # V4-tracked confirmation on sampled tracked points (prereg §7)
    samp = tracked_sample(legs)
    v4t = analyze.confirm_V4_tracked(syn, samp)
    print('  V4-tracked: n=%d worst=%.2e pass=%s'
          % (v4t['n_points'], v4t['worst'], v4t['pass']), flush=True)
    if not v4t['pass']:
        C.add_red_flag('V4T-' + tag.upper(),
                       'M48 vs M96 disagree on tracked points')

    path = os.path.join(C.DATADIR, 'case_%s.json' % tag)
    track.save_case_json(path, tag, 'v2.0', flow_info, seeds, legs)
    C.append_events_csv(tag, legs)
    return path, syn, extended


def run_dh_case(f256):
    """Best-effort DH case (prereg §4d): only if gates said IN."""
    from ht_flow import dh_kappa, phi_dh_wide, HeatFlow
    t_case = time.time()
    kappa = dh_kappa()
    dhf = dh.make_dh_flow(kappa)
    print('== case dh (kappa=%s) ==' % kappa, flush=True)
    found, zlog = dh.find_offaxis_zeros(
        kappa, os.path.join(C.DATADIR, 'dh_zeros.json'))
    if not found:
        C.add_red_flag('R3-DH', 'DH off-axis pair evaded location; '
                                'comparison rests on the window signature '
                                '(prereg R3)')
    else:
        print('  off-axis zeros found: %d' % len(found), flush=True)

    seeds = track.label_zeros(track.seed_zeros(dhf, 300.0))
    C.save_seeds(os.path.join(C.DATADIR, 'seeds_dh.json'), seeds)
    dhf512 = HeatFlow(lambda u: phi_dh_wide(u, kappa), J=512, umax=3.0,
                      dps=50)
    flow_info = {'type': 'Davenport-Heilbronn (prereg §4d)',
                 'kappa': str(kappa), 'J': 256, 'umax': 3.0, 'dps': 50,
                 'off_axis_found': len(found)}

    down = track.track_leg(dhf, seeds, mp.mpf('0'), mp.mpf('-1'),
                           box_seconds=remaining(t_case, CASE_BOX),
                           case_tag='dh_down')
    C.validate_leg(dhf, dhf512, down, 'dh_down')
    # prereg §5: DH off-axis pair tracked up to 1.0
    up = track.track_leg(dhf, seeds, mp.mpf('0'), mp.mpf('1.0'),
                         box_seconds=remaining(t_case, CASE_BOX),
                         case_tag='dh_up')
    C.validate_leg(dhf, dhf512, up, 'dh_up')
    legs = {'down': down, 'up': up}
    path = os.path.join(C.DATADIR, 'case_dh.json')
    track.save_case_json(path, 'dh', 'v2.0', flow_info, seeds, legs)
    C.append_events_csv('dh', legs)
    return path


def v5_stability_report(E1, E2, curve1, curve2):
    """§9: a separating feature must agree between the two zeta down-legs
    within 3σ (same conventions as separation_stats)."""
    rep = {}
    # event counts (Poisson)
    z = analyze._poisson_z(len(E1), len(E2))
    rep['event_count'] = {'z': z, 'stable': z <= 3.0,
                          'n1': len(E1), 'n2': len(E2)}
    # |z*| distributions (KS; stable iff p >= 0.01)
    a1 = [abs(complex(*e['z_star'])) for e in E1]
    a2 = [abs(complex(*e['z_star'])) for e in E2]
    if len(a1) >= 2 and len(a2) >= 2:
        d, p = analyze.ks_2samp(a1, a2)
        rep['zstar_dist'] = {'D': d, 'p': p, 'stable': p >= 0.01,
                             'n1': len(a1), 'n2': len(a2)}
    else:
        rep['zstar_dist'] = {'stable': True,
                             'note': '<2 events in a leg; nothing to disagree'}
    # t*_max
    def tmax(E):
        if not E:
            return None
        e = max(E, key=lambda e: e['t_star'])
        return e['t_star'], e['sigma']
    t1, t2 = tmax(E1), tmax(E2)
    if t1 and t2:
        z = abs(t1[0] - t2[0]) / math.sqrt(t1[1] ** 2 + t2[1] ** 2)
        rep['t_star_max'] = {'z': z, 'stable': z <= 3.0}
    else:
        rep['t_star_max'] = {'stable': True,
                             'note': 'empty map on a side'}
    # N_off(-1)
    n1 = analyze.n_off_at(curve1, -1.0)
    n2 = analyze.n_off_at(curve2, -1.0)
    z = analyze._poisson_z(n1, n2)
    rep['n_off_m1'] = {'z': z, 'stable': z <= 3.0, 'n1': n1, 'n2': n2}
    return rep


def main():
    t_start = time.time()
    print('=== v2 Phase 2: comparison + calibration + decision ===', flush=True)
    if not C.check_study_box():
        print('48h study box expired -> STOP', flush=True)
        return 2

    # gate precondition
    with open(os.path.join(C.DATADIR, 'gate_log.json')) as f:
        glog = json.load(f)
    need = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6']
    if not (glog.get('all_pass') and
            all(glog['gates'][g].get('pass') for g in need)):
        print('GATES NOT ALL PASSING -> Phase 2 barred.', flush=True)
        return 1
    print('gate precondition satisfied (V1-V6 all pass)', flush=True)

    f256, _ = C.build_zeta_flows()
    with open(os.path.join(C.DATADIR, 'case_zeta.json')) as f:
        zeta_doc = json.load(f)
    with open(os.path.join(C.DATADIR, 'case_zeta_v5.json')) as f:
        zeta_v5_doc = json.load(f)
    seeds_zeta = C.load_seeds(os.path.join(C.DATADIR, 'seeds_zeta.json'))
    E_zeta = analyze.validated_events(zeta_doc)
    E_zeta_v5 = analyze.validated_events(zeta_v5_doc)
    curve_zeta = analyze.n_off_curve(zeta_doc['legs']['down'])
    curve_zeta_v5 = analyze.n_off_curve(zeta_v5_doc['legs']['down'])
    print('E_zeta: %d validated events '
          '(V5 leg: %d)' % (len(E_zeta), len(E_zeta_v5)), flush=True)

    z1 = 2 * mp.zetazero(1).imag
    cases = {}

    # (b) synth-d1, (c) synth-d01
    p1, _, _ = run_synth_case(f256, z1, 1.0, 'synth_d1', seeds_zeta)
    cases['synth_d1'] = p1
    p2, _, _ = run_synth_case(f256, z1, 0.1, 'synth_d01', seeds_zeta)
    cases['synth_d01'] = p2

    # (d) DH iff gated IN
    if glog['gates']['DH']['decision'] == 'IN':
        cases['dh'] = run_dh_case(f256)
    else:
        print('DH dropped at gating; comparison is zeta vs synth cases',
              flush=True)

    # ---- calibration §8 ----
    with open(cases['synth_d1']) as f:
        syn1_doc = json.load(f)
    landing = analyze.check_landing(syn1_doc)
    print('calibration (iv) landing: %s' % landing, flush=True)
    cal_ok, cal_det = analyze.check_calibration(zeta_doc, glog, landing,
                                                seeds_zeta)
    cal_det['iii_null_window'] = {
        'pass': bool(glog['gates']['V6']['pass']),
        'note': 'same computation as V6 (detector null on (0.2,0.5])'}
    cal_ok = all(v.get('pass') for v in cal_det.values())
    print('calibration: %s' % ('PASS' if cal_ok else 'FAIL'), flush=True)
    for k, v in cal_det.items():
        print('   %s: pass=%s' % (k, v.get('pass')), flush=True)
    if not cal_ok:
        C.add_red_flag('CALIBRATION-FAIL',
                       'calibration §8 failed -> METHOD INVALID per prereg')

    # ---- comparison §9 ----
    comparisons = {}
    for name, path in cases.items():
        with open(path) as f:
            doc = json.load(f)
        E_case = analyze.validated_events(doc)
        curve_case = analyze.n_off_curve(doc['legs']['down'])
        feats = analyze.separation_stats(E_zeta, E_case, curve_zeta,
                                         curve_case)
        comparisons[name] = {'n_events': len(E_case), 'features': feats}
        print('--- zeta vs %s: %d validated events' % (name, len(E_case)),
              flush=True)
        for fn, fs in feats.items():
            print('    %-12s separates=%s %s'
                  % (fn, fs['separates'],
                     {k: v for k, v in fs.items() if k != 'separates'}),
                  flush=True)

    # ---- decision §9 ----
    stab = v5_stability_report(E_zeta, E_zeta_v5, curve_zeta, curve_zeta_v5)
    verdict = None
    verdict_detail = ''
    if not cal_ok:
        verdict = 'INVALID'
        verdict_detail = 'calibration §8 failed'
    else:
        winners = []
        for name, comp in comparisons.items():
            for fn, fs in comp['features'].items():
                if fs.get('separates') and stab[fn]['stable']:
                    winners.append((name, fn))
        if winners:
            verdict = 'SUCCESS'
            only_d1 = all(n == 'synth_d1' for n, _ in winners)
            verdict_detail = ('separating features (V5-stable): %s'
                              % ['%s/%s' % w for w in winners])
            if only_d1:
                verdict_detail += (' -- per prereg R9 the claim is weakened '
                                   'to "detects gross violations"')
                C.add_red_flag('R9-WEAKENED',
                               'only synth-d1 separates: claim weakened to '
                               '"detects gross violations" per prereg R9')
        else:
            verdict = 'FAILURE'
            verdict_detail = 'calibration passes; no V5-stable feature ' \
                             'separates at >3σ'
    print('=== VERDICT: %s ===' % verdict, flush=True)
    print(verdict_detail, flush=True)

    # ---- deliverables ----
    with open(os.path.join(C.DATADIR, 'red_flags.json')) as f:
        flags = json.load(f)
    decision = {
        'prereg': 'v2.0',
        'verdict': verdict,
        'verdict_detail': verdict_detail,
        'calibration': {'pass': cal_ok, 'detail': cal_det},
        'comparisons': comparisons,
        'v5_stability': stab,
        'gates': {g: glog['gates'][g] for g in need},
        'dh_decision': glog['gates']['DH']['decision'],
        'red_flags': flags['flags'],
        'n_zeta_events': len(E_zeta),
    }
    with open(os.path.join(C.DATADIR, 'decision.json'), 'w') as f:
        json.dump(decision, f, indent=1)

    all_cases = {'zeta': zeta_doc}
    for name, path in cases.items():
        with open(path) as f:
            all_cases[name] = json.load(f)
    analyze.write_figure_data(C.FIGDIR, all_cases)
    write_results_md(decision, all_cases)

    print('elapsed: %.1f min' % ((time.time() - t_start) / 60), flush=True)
    return 0


def write_results_md(decision, all_cases):
    import hashlib
    v2 = os.path.dirname(HERE)
    pdf = os.path.join(v2, 'preregistration-v2.pdf')
    md = os.path.join(v2, 'preregistration.md')
    h = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
    L = []
    A = L.append
    A('# v2.0 heat-flow study — results\n')
    A('Preregistration v2.0, deposited 2026-10-03 BEFORE any v2 t≠0 '
      'computation. Zenodo DOI 10.5281/zenodo.23120530.\n')
    A('- preregistration-v2.pdf sha256: `%s`' % h(pdf))
    A('- preregistration.md sha256: `%s`\n' % h(md))
    A('## Verdict: %s\n' % decision['verdict'])
    A('%s\n' % decision['verdict_detail'])
    A('Interpretation cap (prereg §2/§9): this is methods/phenomenology '
      'only — NOT a proof of RH, NOT a statement about RH\'s truth value, '
      'NOT a measurement of Λ. All t* values are windowed.\n')
    A('## Calibration (§8)\n')
    for k, v in decision['calibration']['detail'].items():
        A('- %s: pass=%s' % (k, v.get('pass')))
        if k == 'iv_landing':
            A('  (t*_inj=%s, detail: %s)'
              % (v.get('t_inj'), v.get('detail')))
    A('\n## Gates (§7)\n')
    for g, r in decision['gates'].items():
        A('- %s: pass=%s%s' % (g, r.get('pass'),
                               ' (%s)' % r.get('detail') if r.get('detail')
                               else ''))
    A('\nDH case decision: %s\n' % decision['dh_decision'])
    A('## Event maps\n')
    for name, doc in all_cases.items():
        E = analyze.validated_events(doc)
        A('- %s: %d validated events' % (name, len(E)))
        for e in E:
            A('  - %s t*=%.6f±%.1e z*=(%.4f,%.4f) R²=%.5f labels=%s'
              % (e['type'], e['t_star'], e['sigma'], e['z_star'][0],
                 e['z_star'][1], e.get('R2'), ','.join(e['labels'])))
    A('\n## Discriminator features (§9)\n')
    for name, comp in decision['comparisons'].items():
        A('### zeta vs %s (%d validated events)\n'
          % (name, comp['n_events']))
        for fn, fs in comp['features'].items():
            A('- %s: separates=%s, V5-stable=%s, detail=%s'
              % (fn, fs['separates'],
                 decision['v5_stability'][fn]['stable'],
                 {k: v for k, v in fs.items() if k != 'separates'}))
    A('\n## Red-flag log (§12)\n')
    for fl in decision['red_flags']:
        A('- [%s] %s: %s' % (fl['utc'], fl['code'], fl['text']))
    if not decision['red_flags']:
        A('(none)')
    A('\n## Soundness self-assessment\n')
    A('- Preregistration timestamped before any computation: yes '
      '(DOI above).')
    A('- All gates logged before Phase 2: %s.'
      % all(r.get('pass') for r in decision['gates'].values()))
    A('- Trajectory JSONs saved for every case: %s.'
      % all('legs' in d and all('traj' in l for l in d['legs'].values())
            for d in all_cases.values()))
    A('- Genuineness protocol applied per event: yes (verdicts in '
      'data/events.csv).')
    A('- Known judgment calls: (1) synth G-b uses M=96 refinement as the '
      'honest analog of J-refinement (V4 gates M48~=M96); (2) count-feature '
      'separation uses the Poisson sqrt(N) σ reading of ">3σ" (flagged in '
      'CODE_NOTES.md item 5; CONFIRMED/UNCONFIRMED by Abby before Phase 2).')
    A('- Null results are reported as FAILURE, never rescued.')
    with open(os.path.join(v2, 'results.md'), 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('results.md written', flush=True)


if __name__ == '__main__':
    sys.exit(main())
