#!/usr/bin/env python3
"""run_phase1.py — prereg v2.0 §9 Phase 1: zeta characterization + V5.

Produces the validated event map E_zeta and artifact tally A_zeta
(descriptive; no discriminator claim). Then V5: repeat the down-leg with
seeds perturbed by 1e-6 and require the validated event maps to agree
(same events within 3σ, none gained/lost). The V5 result is appended to
data/gate_log.json; Phase 2 is barred unless every gate passes.
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

ZETA_DOWN_BOX = 6 * 3600  # prereg §10


def main():
    t_start = time.time()
    print('=== v2 Phase 1: zeta characterization ===', flush=True)
    if not C.check_study_box():
        print('48h study box expired -> STOP', flush=True)
        return 2
    f256, f512 = C.build_zeta_flows()
    seeds = C.load_seeds(os.path.join(C.DATADIR, 'seeds_zeta.json'))
    print('loaded %d zeta seeds' % len(seeds), flush=True)

    print('-- zeta down-leg [0, -1] (6h box) --', flush=True)
    down = track.track_leg(f256, seeds, mp.mpf('0'), mp.mpf('-1'),
                           box_seconds=ZETA_DOWN_BOX, case_tag='zeta_down')
    print('   reached t=%.4f, steps=%d, box_expired=%s, track_fails=%d'
          % (down['t_reached'], down['n_steps'], down['box_expired'],
             len(down['track_fails'])), flush=True)
    if down['box_expired']:
        C.add_red_flag('BOX-ZETA-DOWN',
                       'zeta down-leg hit the 6h box; partial map reported '
                       'per prereg §10 (scope cut, not a failure)')
    C.validate_leg(f256, f512, down, 'zeta_down')

    print('-- zeta up-leg [0, 0.5] --', flush=True)
    up = track.track_leg(f256, seeds, mp.mpf('0'), mp.mpf('0.5'),
                         case_tag='zeta_up')
    C.validate_leg(f256, f512, up, 'zeta_up')

    legs = {'down': down, 'up': up}
    track.save_case_json(os.path.join(C.DATADIR, 'case_zeta.json'), 'zeta',
                         'v2.0', C.zeta_flow_info(), seeds, legs)
    print('case -> data/case_zeta.json', flush=True)
    C.append_events_csv('zeta', legs)

    # artifact tally (prereg §6: rate > 0.5 -> INVALID apparatus)
    evs = [e for leg in legs.values() for e in leg['events']]
    tally = analyze.artifact_tally([track.event_to_jsonable(e)
                                    for e in evs])
    print('artifact tally: %s, rate = %.3f'
          % (tally['tally'], tally['artifact_rate']), flush=True)
    if tally['artifact_rate'] > 0.5:
        C.add_red_flag('ARTIFACT-RATE',
                       'zeta artifact rate %.3f > 0.5 -> INVALID (apparatus) '
                       'per prereg §6' % tally['artifact_rate'])
    with open(os.path.join(C.DATADIR, 'phase1_summary.json'), 'w') as f:
        json.dump({'artifact_tally': tally,
                   'down': {'t_reached': down['t_reached'],
                            'n_steps': down['n_steps'],
                            'box_expired': down['box_expired'],
                            'track_fails': down['track_fails']},
                   'n_events_down': len(down['events']),
                   'n_events_up': len(up['events'])}, f, indent=1)

    # ---- V5: seed independence (prereg §7) ----
    print('-- V5: perturbed-seed repeat of the down-leg --', flush=True)
    pert = {lab: z + mp.mpc('1e-6', '1e-6') for lab, z in seeds.items()}
    down2 = track.track_leg(f256, pert, mp.mpf('0'), mp.mpf('-1'),
                            box_seconds=ZETA_DOWN_BOX, case_tag='zeta_down_V5')
    C.validate_leg(f256, f512, down2, 'zeta_down_V5')
    track.save_case_json(os.path.join(C.DATADIR, 'case_zeta_v5.json'),
                         'zeta-V5', 'v2.0', C.zeta_flow_info(), pert,
                         {'down': down2})

    doc_a = {'legs': {'down': {'events':
                               [track.event_to_jsonable(e)
                                for e in down['events']]}}}
    doc_b = {'legs': {'down': {'events':
                               [track.event_to_jsonable(e)
                                for e in down2['events']]}}}
    E_a = analyze.validated_events(doc_a)
    E_b = analyze.validated_events(doc_b)
    cmp_res = analyze.compare_event_maps(E_a, E_b)
    v5 = {'gate': 'V5', 'pass': bool(cmp_res['pass']),
          'detail': cmp_res['detail'],
          'n_events_leg1': len(E_a), 'n_events_leg2': len(E_b)}
    print('   V5: %s (%s)' % (v5['pass'], v5['detail']), flush=True)

    # append V5 to the gate log and recompute all_pass
    gp = os.path.join(C.DATADIR, 'gate_log.json')
    with open(gp) as f:
        glog = json.load(f)
    glog['gates']['V5'] = v5
    glog['all_pass'] = all(r.get('pass') for r in glog['gates'].values())
    with open(gp, 'w') as f:
        json.dump(glog, f, indent=1)
    print('gate log updated; all_pass = %s' % glog['all_pass'], flush=True)

    # figure data for zeta so far
    with open(os.path.join(C.DATADIR, 'case_zeta.json')) as f:
        zdoc = json.load(f)
    analyze.write_figure_data(C.FIGDIR, {'zeta': zdoc})

    print('elapsed: %.1f min' % ((time.time() - t_start) / 60), flush=True)
    if not v5['pass']:
        C.add_red_flag('V5-FAIL', 'seed-independence failed: ' +
                       cmp_res['detail'])
        print('V5 FAILED -> STOP. Phase 2 barred.', flush=True)
        return 1
    print('Phase 1 complete. E_zeta validated events: %d. '
          'Proceed to run_phase2.py.' % len(E_a), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
