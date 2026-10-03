#!/usr/bin/env python3
"""Shared helpers for the v2 heat-flow runners (preregistration v2.1)."""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
V2DIR = os.path.dirname(HERE)
DATADIR = os.path.join(V2DIR, 'data')
FIGDIR = os.path.join(V2DIR, 'figures')
os.makedirs(DATADIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

from mpmath import mp  # noqa: E402
import track  # noqa: E402
import analyze  # noqa: E402

STUDY_BOX_HOURS = 48.0  # prereg §10: total study box


def build_zeta_flows():
    from ht_flow import HeatFlow, phi_zeta
    # v2.1 repair: dps 50->80 (rounding killed Newton in flat bands),
    # umax 3.0->1.2 (Phi mass is in [0,1]; denser nodes fix oscillatory
    # truncation: J256->512 diff 2e-71 uniform at z<=300, was 4e-22).
    f256 = HeatFlow(phi_zeta, J=256, umax=1.2, dps=80)
    f512 = HeatFlow(phi_zeta, J=512, umax=1.2, dps=80)
    return f256, f512


def zeta_flow_info():
    return {'phi': 'phi_zeta (prereg §4a)', 'J': 256, 'umax': 1.2, 'dps': 80,
            'note': 'v2.1 repaired instrument (see preregistration-v2.1 §2)'}


def save_seeds(path, seeds):
    doc = {lab: [float(mp.re(z)), float(mp.im(z))]
           for lab, z in seeds.items()}
    with open(path, 'w') as f:
        json.dump(doc, f)


def load_seeds(path):
    with open(path) as f:
        doc = json.load(f)
    return {lab: mp.mpc(v[0], v[1]) for lab, v in doc.items()}


def stamp_study_start():
    """Record the first-v2-t!=0-evaluation time (prereg §10: 48h box)."""
    p = os.path.join(DATADIR, 'study_state.json')
    if not os.path.exists(p):
        with open(p, 'w') as f:
            json.dump({'first_tnonzero_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                                          time.gmtime()),
                       'prereg': 'v2.1'}, f)


def check_study_box():
    p = os.path.join(DATADIR, 'study_state.json')
    if not os.path.exists(p):
        return True
    with open(p) as f:
        st = json.load(f)
    t0 = time.mktime(time.strptime(st['first_tnonzero_utc'],
                                   '%Y-%m-%dT%H:%M:%SZ'))
    # strptime gives local; study_state was written with gmtime -> use calendar
    import calendar
    t0 = calendar.timegm(time.strptime(st['first_tnonzero_utc'],
                                        '%Y-%m-%dT%H:%M:%SZ'))
    elapsed_h = (time.time() - t0) / 3600.0
    print('study box: %.1f / %.0f h elapsed' % (elapsed_h, STUDY_BOX_HOURS),
          flush=True)
    return elapsed_h < STUDY_BOX_HOURS


def validate_leg(flow_a, flow_b, leg, tag):
    """Run the §6 genuineness protocol on every event of a raw leg
    (in place). flow_a/flow_b: (256, 512) for quadrature cases,
    (M48, M96) for synthetic cases (M-refinement is the honest analog of
    J-refinement; V4 gates exactly this equivalence)."""
    for ev in leg['events']:
        track.validate_event(flow_a, flow_b, ev, leg['traj'])
    n = {'genuine': 0, 'suspect': 0, 'artifact': 0, 'logged': 0}
    for ev in leg['events']:
        n[ev.get('verdict', '?')] = n.get(ev.get('verdict', '?'), 0) + 1
    print('  [%s] events: %d -> %s' % (tag, len(leg['events']), n), flush=True)


def append_events_csv(case, legs):
    """Append validated, serialized events to data/events.csv."""
    import csv
    path = os.path.join(DATADIR, 'events.csv')
    new = not os.path.exists(path)
    cols = ['case', 'leg', 'type', 't_detect', 't_star', 'sigma', 'R2', 'C',
            'z_star_re', 'z_star_im', 'n_fit_pts', 'labels',
            'G_a', 'G_b', 'G_c', 'verdict', 'artifact_reason']
    with open(path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        if new:
            w.writeheader()
        for leg_name, leg in legs.items():
            for ev in leg['events']:
                row = analyze.event_row(case, leg_name,
                                        track.event_to_jsonable(ev))
                w.writerow({c: row.get(c, '') for c in cols})


def red_flags():
    """Load-or-create the study red-flag log (prereg §12, results §11)."""
    p = os.path.join(DATADIR, 'red_flags.json')
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f), p
    log = {'prereg': 'v2.1', 'flags': []}
    return log, p


def add_red_flag(code, text):
    log, p = red_flags()
    log['flags'].append({'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                              time.gmtime()),
                         'code': code, 'text': text})
    with open(p, 'w') as f:
        json.dump(log, f, indent=1)
    print('RED FLAG [%s]: %s' % (code, text), flush=True)
