#!/usr/bin/env python3
"""
analyze.py — verification gates V1-V6 (§7), calibration checks (§8),
event-map assembly, discriminator statistics (§9), artifact tally,
events CSV and figure data (§11) for preregistration v2.0.

Conventions (prereg literal):
  - N_off counts zeros with |Im z| > 1e-6.
  - KS test on |z*|: two-sample, alpha = 0.01; <5 events in either map =>
    non-separating (never a success).
  - >3sigma separation, implemented as:
      |E| counts : Poisson z = |E1-E2| / sqrt(E1+E2) > 3
      t*_max     : |t1-t2| / sqrt(s1^2+s2^2) > 3   (fit sigmas)
      N_off(-1)  : Poisson z = |n1-n2| / sqrt(n1+n2) > 3
      |z*| dist  : KS p < 0.01 (not a sigma rule; prereg specifies KS)
    NOTE (underspecification, flagged): the preregistration does not define
    the sigma for count features; the Poisson-sqrt(N) choice above is the
    natural reading and must be confirmed before Phase 2 runs.
  - V5 stability: a separating feature must agree between the two zeta
    down-legs within 3sigma under the same conventions.
"""
import csv
import json
import math
import os
import time
from mpmath import mp

import track

OFF_TOL = mp.mpf('1e-6')
KS_ALPHA = 0.01
KS_MIN_EVENTS = 5
SEPARATION_SIGMAS = 3.0


# ---------- loading ----------

def load_case_json(path):
    with open(path) as f:
        return json.load(f)


def validated_events(case_doc):
    """All events with verdict == 'genuine' across legs."""
    out = []
    for leg in case_doc['legs'].values():
        for ev in leg['events']:
            if ev.get('verdict') == 'genuine':
                out.append(ev)
    return out


def all_triggers(case_doc):
    out = []
    for leg in case_doc['legs'].values():
        out.extend(leg['events'])
    return out


# ---------- two-sample Kolmogorov-Smirnov ----------

def ks_2samp(x, y):
    """Returns (D, p) with the asymptotic Kolmogorov distribution
    (finite-sample correction). Requires len >= 1 each; the prereg's
    >= 5 rule is enforced by the caller."""
    nx, ny = len(x), len(y)
    xs = sorted(x)
    ys = sorted(y)
    i = j = 0
    d = 0.0
    # merge walk over sorted unique values
    vals = sorted(set(xs) | set(ys))
    for v in vals:
        while i < nx and xs[i] <= v:
            i += 1
        while j < ny and ys[j] <= v:
            j += 1
        d = max(d, abs(i / nx - j / ny))
    en = nx * ny / (nx + ny)
    lam = (math.sqrt(en) + 0.12 + 0.11 / math.sqrt(en)) * d
    # Kolmogorov survival function
    p = 0.0
    k = 1
    while True:
        term = (-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam)
        p += term
        if abs(term) < 1e-12 or k > 100:
            break
        k += 1
    p = max(0.0, min(1.0, 2 * p))
    return d, p


# ---------- N_off(t) ----------

def n_off_curve(leg):
    """[(t, N_off)] with N_off = #{labels: |Im z| > 1e-6} per step."""
    traj = leg['traj']
    # union of step times
    ts = sorted({p[0] for pts in traj.values() for p in pts})
    # index per label
    curve = []
    for t in ts:
        n = 0
        for pts in traj.values():
            # nearest recorded point at-or-before t
            best = None
            for p in pts:
                if p[0] <= t + 1e-12:
                    best = p
                else:
                    break
            if best is not None and abs(best[2]) > float(OFF_TOL):
                n += 1
        curve.append((t, n))
    return curve


def n_off_at(curve, t):
    return min(curve, key=lambda tn: abs(tn[0] - t))[1]


# ---------- artifact tally ----------

def artifact_tally(events):
    tally = {'genuine': 0, 'suspect': 0, 'artifact': 0, 'logged': 0}
    reasons = []
    for ev in events:
        v = ev.get('verdict', 'unvalidated')
        tally[v] = tally.get(v, 0) + 1
        if v in ('suspect', 'artifact') and ev.get('artifact_reason'):
            reasons.append({'type': ev.get('type'), 'labels': ev.get('labels'),
                            't_detect': ev.get('t_detect'),
                            'reason': ev['artifact_reason']})
    n = tally['genuine'] + tally['artifact']
    rate = tally['artifact'] / n if n else 0.0
    return {'tally': tally, 'artifact_rate': rate, 'reasons': reasons}


# ---------- gates V1-V6 (prereg §7) ----------

def _test_points_zeta():
    # 10 (z,t); 4 with |z| in [180,280] (V1 requirement)
    zs = [10.0, 37.0, 80.0, 150.0, 190.0, 210.0, 240.0, 275.0, 60.0, 120.0]
    ts = [-0.5, -0.1, 0.0, 0.2, 0.4, -0.5, -0.1, 0.0, 0.2, 0.4]
    return list(zip(zs, ts))


def gate_V1(flow256, flow512):
    """J=256 -> J=512 changes H_t(z) by < 1e-60 at 10 (z,t), incl. 4 with
    |z| in [180,280].

    v2.1 recalibration: v2.0's 1e-35 was miscalibrated for the wide window
    (measured 3.79e-22 at z=275, umax=3.0, dps=50). With the repaired
    instrument (umax=1.2, dps=80) the measured J-refinement is ~2e-71
    uniform in z; threshold 1e-60 keeps 11 orders of margin.
    """
    worst = 0.0
    det = []
    for zr, tr in _test_points_zeta():
        d = abs(flow256.H(zr, tr) - flow512.H(zr, tr))
        worst = max(worst, float(d))
        det.append({'z': zr, 't': tr, 'dH': float(d)})
    return {'gate': 'V1', 'worst': worst, 'threshold': 1e-60,
            'pass': worst < 1e-60, 'points': det}


def gate_V2(flow):
    """|H_0(2γ_n)| < 1e-30 for n=1..5; 14-digit agreement with 2γ_n."""
    det = []
    worst_h = 0.0
    worst_rel = 0.0
    for n in range(1, 6):
        g = mp.zetazero(n).imag
        z0 = 2 * g
        h = abs(flow.H(z0, 0))
        worst_h = max(worst_h, float(h))
        # polish and compare
        z, _ = track.newton(flow, z0, mp.mpf('0'))
        rel = float(abs(z - z0) / abs(z0)) if z is not None else float('inf')
        worst_rel = max(worst_rel, rel)
        det.append({'n': n, 'H0': float(h), 'rel_err': rel})
    return {'gate': 'V2', 'worst_H0': worst_h, 'threshold_H0': 1e-30,
            'worst_rel': worst_rel, 'threshold_rel': 1e-14,
            'pass': worst_h < 1e-30 and worst_rel < 1e-14, 'points': det}


def gate_V2b(flow, zmax=300.0):
    """Seeding completeness (v2.1, new): the grid-seeded Newton finder
    recovers every known zeta zero with |2γ_n| <= zmax.

    v2.0 failed here: at dps=50, the quadrature suffered ~35-digit
    cancellation in z-dependent high-z bands, seeding found only 60 of
    ~104 expected zeros, and Newton converged to displaced artifact roots
    (up to ~5 units from any true zero). v2.1 requires: count matches the
    reference count exactly, every seeded zero within 0.01 of a reference
    ±2γ_n, and worst |H_0| < 1e-25.
    """
    import track as _track
    zeros = _track.seed_zeros(flow, zmax=zmax)
    # reference zeros: ±2γ_n with |2γ_n| <= zmax
    refs = []
    n = 1
    while True:
        g = float(mp.zetazero(n).imag)
        z0 = 2 * g
        if z0 > zmax:
            break
        refs.append(z0)
        refs.append(-z0)
        n += 1
    n_ref = len(refs)
    # nearest-reference distance for each seeded zero
    worst_disp = 0.0
    worst_h = 0.0
    for z in zeros:
        zf = float(mp.re(z))
        worst_h = max(worst_h, float(abs(flow.H(z, mp.mpf('0')))))
        dmin = min(abs(zf - r) for r in refs) if refs else float('inf')
        worst_disp = max(worst_disp, dmin)
    # reference coverage: every ref within 0.01 of a seeded zero
    zf_list = [float(mp.re(z)) for z in zeros]
    worst_miss = 0.0
    for r in refs:
        dmin = min(abs(zf - r) for zf in zf_list) if zf_list else float('inf')
        worst_miss = max(worst_miss, dmin)
    ok = (len(zeros) == n_ref and worst_disp < 0.01 and worst_miss < 0.01
          and worst_h < 1e-25)
    return {'gate': 'V2b', 'n_seeded': len(zeros), 'n_reference': n_ref,
            'worst_displacement': worst_disp,
            'worst_reference_miss': worst_miss, 'worst_H0': worst_h,
            'thresholds': {'displacement': 0.01, 'miss': 0.01, 'H0': 1e-25},
            'pass': bool(ok)}


def gate_V3(flow, case_tag):
    """|dH/dt + d2H/dz2| < 1e-60 at 10 (z,t), analytic dH/dt.

    v2.1 repair: v2.0 used a central finite difference whose O(dt^2)
    truncation (1.25e-19 at dt=1e-6) swamped the true residual (~1e-42
    estimated). pde_residual_exact uses the analytic time derivative, so
    the gate measures the true quadrature PDE residual (~1e-85 at dps=80).
    """
    worst = 0.0
    det = []
    for zr, tr in _test_points_zeta():
        r = flow.pde_residual_exact(mp.mpf(str(zr)), mp.mpf(str(tr)))
        worst = max(worst, float(r))
        det.append({'z': zr, 't': tr, 'res': float(r)})
    return {'gate': 'V3', 'case': case_tag, 'worst': worst,
            'threshold': 1e-60, 'pass': worst < 1e-60, 'points': det}


def gate_V4(syn48, w):
    """M=48 vs M=96 agree < 1e-25 at 10 test (z,t); H_0^syn(w)=0 < 1e-25;
    PDE residual (analytic dH/dt) < 1e-60. (Post-tracking confirmatory
    re-run on tracked points: see confirm_V4_tracked.)

    v2.1: PDE via pde_residual_exact; v2.0's finite-difference version
    measured its own truncation (marginal 1.01e-20 'fail')."""
    from synthetic import SyntheticFlow
    syn96 = SyntheticFlow(syn48.base, syn48.z1, syn48.delta, M=96)
    worst = 0.0
    det = []
    for zr, tr in _test_points_zeta():
        d = abs(syn48.H(zr, tr) - syn96.H(zr, tr))
        worst = max(worst, float(d))
        det.append({'z': zr, 't': tr, 'dH': float(d)})
    hc = syn48.construction_check()
    pde = syn48.pde_residual_exact(mp.mpf('37'), mp.mpf('0.1'))
    return {'gate': 'V4', 'worst_Mdiff': worst, 'threshold_Mdiff': 1e-25,
            'construction': float(hc), 'threshold_construction': 1e-25,
            'pde': float(pde), 'threshold_pde': 1e-60,
            'pass': worst < 1e-25 and float(hc) < 1e-25 and float(pde) < 1e-60,
            'points': det}


def confirm_V4_tracked(syn48, traj_sample):
    """Post-tracking confirmatory: M48 vs M96 < 1e-25 at sampled tracked
    (z,t) points. traj_sample: list of (z_re, z_im, t)."""
    from synthetic import SyntheticFlow
    syn96 = SyntheticFlow(syn48.base, syn48.z1, syn48.delta, M=96)
    worst = 0.0
    for zr, zi, tr in traj_sample:
        z = mp.mpc(zr, zi)
        d = abs(syn48.H(z, tr) - syn96.H(z, tr))
        worst = max(worst, float(d))
    return {'check': 'V4-tracked', 'n_points': len(traj_sample),
            'worst': worst, 'threshold': 1e-25, 'pass': worst < 1e-25}


def compare_event_maps(E_a, E_b):
    """V5 comparison: same validated events within 3σ, none gained/lost."""
    unmatched_b = list(E_b)
    pairs = []
    for ea in E_a:
        best, best_key = None, None
        for eb in unmatched_b:
            if eb['type'] != ea['type']:
                continue
            dt = abs(eb['t_star'] - ea['t_star'])
            if dt >= 3 * math.sqrt(ea['sigma'] ** 2 + eb['sigma'] ** 2):
                continue
            za = complex(*ea['z_star'])
            zb = complex(*eb['z_star'])
            if abs(za - zb) >= 0.05:
                continue
            key = (dt, abs(za - zb))
            if best is None or key < best_key:
                best, best_key = eb, key
        if best is None:
            return {'pass': False,
                    'detail': 'event in A unmatched: %s t*=%.4f' %
                              (ea['type'], ea['t_star'])}
        pairs.append((ea, best))
        unmatched_b.remove(best)
    if unmatched_b:
        return {'pass': False,
                'detail': '%d event(s) in B unmatched' % len(unmatched_b)}
    return {'pass': True, 'detail': '%d events matched' % len(pairs)}


def gate_V6(flow, seeds):
    """Full event detector on ζ over (0.2, 0.5]; zero events required.
    Seeds are tracked 0 -> 0.2 first, then the detector leg 0.2 -> 0.5."""
    pre = track.track_leg(flow, seeds, mp.mpf('0'), mp.mpf('0.2'),
                          case_tag='zeta_V6pre')
    seeds02 = {lab: pts[-1][1] for lab, pts in pre['traj'].items()}
    leg = track.track_leg(flow, seeds02, mp.mpf('0.2'), mp.mpf('0.5'),
                          case_tag='zeta_V6')
    n_trig = len(leg['events'])
    return {'gate': 'V6', 'n_triggers': n_trig, 'pass': n_trig == 0,
            'detail': 'detector triggers in (0.2,0.5]: %d (require 0)' %
                      n_trig,
            'triggers': [track.event_to_jsonable(e) for e in leg['events']]}


def log_gates(path, results):
    doc = {'prereg': 'v2.0', 'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'gates': results,
           'all_pass': all(r.get('pass') for r in results.values())}
    with open(path, 'w') as f:
        json.dump(doc, f, indent=1)
    return doc['all_pass']


# ---------- calibration §8 ----------

def check_calibration(case_zeta, gate_log, syn1_landing, seeds_zeta):
    """Returns (pass, detail). syn1_landing: dict from check_landing()."""
    det = {}
    # (i) N_off(0) == 0 in |z| <= 300
    n_off0 = sum(1 for v in seeds_zeta.values()
                 if abs(mp.im(mp.mpc(v[0], v[1]))) > float(OFF_TOL))
    det['i_N_off_0'] = {'value': n_off0, 'pass': n_off0 == 0}
    # (ii) gates logged and all pass
    det['ii_gates'] = {'pass': bool(gate_log.get('all_pass'))}
    # (iii) zero validated events in (0.2, 0.5] — same computation as V6
    det['iii_null_window'] = {'pass': True,  # set by runner from V6
                              'note': 'same computation as V6'}
    # (iv) apparatus validation: synth-δ1 landing
    det['iv_landing'] = syn1_landing
    ok = all(v.get('pass') for v in det.values())
    return ok, det


def check_landing(case_syn1, injected_labels=('w1', 'w2', 'w3', 'w4'),
                  t_max=2.0):
    """Calibration §8(iv): the injected quadruplet lands via Type-C at
    0 < t*_inj <= 2.0, im²-fit R² >= 0.999, N_off drops by exactly 4."""
    evs = [e for e in validated_events(case_syn1)
           if e['type'] == 'C'
           and all(lab in injected_labels for lab in e['labels'])]
    if not evs:
        return {'pass': False,
                'detail': 'no validated Type-C landing of injected labels'}
    ts = [e['t_star'] for e in evs]
    sigs = [e['sigma'] for e in evs]
    # inverse-variance weighted mean
    wts = [1 / (s ** 2) for s in sigs]
    t_inj = sum(t * w for t, w in zip(ts, wts)) / sum(wts)
    sig_inj = math.sqrt(1 / sum(wts))
    spread_ok = all(abs(t - t_inj) < 3 * math.sqrt(s ** 2 + sig_inj ** 2)
                    for t, s in zip(ts, sigs))
    r2_ok = all(e['R2'] >= 0.999 for e in evs)
    in_window = 0 < t_inj <= t_max
    # N_off drop across the landing: compare curve just before/after
    leg = None
    for lname, l in case_syn1['legs'].items():
        ev_here = [e for e in l['events']
                   if e.get('verdict') == 'genuine' and e['type'] == 'C'
                   and all(x in injected_labels for x in e['labels'])]
        if ev_here:
            leg = l
            break
    drop = None
    if leg is not None:
        curve = n_off_curve(leg)
        n_lo = n_off_at(curve, t_inj - 5 * sig_inj)
        n_hi = n_off_at(curve, t_inj + 5 * sig_inj)
        drop = n_lo - n_hi
    drop_ok = (drop == 4)
    ok = spread_ok and r2_ok and in_window and drop_ok
    return {'pass': bool(ok), 't_inj': t_inj, 'sigma_inj': sig_inj,
            'n_landings': len(evs), 'spread_ok': spread_ok, 'r2_ok': r2_ok,
            'in_window': in_window, 'n_off_drop': drop, 'drop_ok': drop_ok,
            'detail': 't*_inj=%.6f±%.1e, R2 ok=%s, N_off drop=%s (need 4)' %
                      (t_inj, sig_inj, r2_ok, drop)}


# ---------- discriminator statistics §9 ----------

def _poisson_z(a, b):
    s = a + b
    if s <= 0:
        return 0.0
    return abs(a - b) / math.sqrt(s)


def separation_stats(E_zeta, E_case, curve_zeta_down, curve_case_down):
    """Preregistered features; returns {feature: {separates: bool, ...}}."""
    out = {}
    # |E| event counts (Poisson)
    ez, ec = len(E_zeta), len(E_case)
    z = _poisson_z(ez, ec)
    out['event_count'] = {'zeta': ez, 'case': ec, 'z': z,
                          'separates': z > SEPARATION_SIGMAS}
    # |z*| distribution (KS, alpha=0.01; >=5 events each)
    az = [abs(complex(*e['z_star'])) for e in E_zeta]
    ac = [abs(complex(*e['z_star'])) for e in E_case]
    if len(az) >= KS_MIN_EVENTS and len(ac) >= KS_MIN_EVENTS:
        d, p = ks_2samp(az, ac)
        out['zstar_dist'] = {'D': d, 'p': p, 'separates': p < KS_ALPHA,
                             'n_zeta': len(az), 'n_case': len(ac)}
    else:
        out['zstar_dist'] = {'separates': False,
                             'note': 'insufficient events (<5 in a map); '
                                     'non-separating per prereg §9'}
    # t*_max (latest event time), fit-sigma weighted
    def tmax(E):
        if not E:
            return None
        e = max(E, key=lambda e: e['t_star'])
        return e['t_star'], e['sigma']
    tz, tc = tmax(E_zeta), tmax(E_case)
    if tz and tc:
        zt = abs(tz[0] - tc[0]) / math.sqrt(tz[1] ** 2 + tc[1] ** 2)
        out['t_star_max'] = {'zeta': tz[0], 'case': tc[0], 'z': zt,
                             'separates': zt > SEPARATION_SIGMAS}
    else:
        out['t_star_max'] = {'separates': False,
                             'note': 'empty event map on a side'}
    # N_off(-1)
    nz = n_off_at(curve_zeta_down, -1.0)
    nc = n_off_at(curve_case_down, -1.0)
    zn = _poisson_z(nz, nc)
    out['n_off_m1'] = {'zeta': nz, 'case': nc, 'z': zn,
                       'separates': zn > SEPARATION_SIGMAS}
    return out


def v5_stability(feature_name, val_a, val_b, sigma_a, sigma_b):
    """A separating feature must agree between the two zeta down-legs
    within 3sigma (prereg §9)."""
    z = abs(val_a - val_b) / math.sqrt(sigma_a ** 2 + sigma_b ** 2)
    return {'stable': z <= SEPARATION_SIGMAS, 'z': z}


# ---------- outputs §11 ----------

def write_events_csv(path, rows):
    cols = ['case', 'leg', 'type', 't_detect', 't_star', 'sigma', 'R2', 'C',
            'z_star_re', 'z_star_im', 'n_fit_pts', 'labels',
            'G_a', 'G_b', 'G_c', 'verdict', 'artifact_reason']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, '') for c in cols})


def event_row(case, leg, ev):
    g = ev.get('genuineness', {})
    zs = ev.get('z_star') or [None, None]
    return {
        'case': case, 'leg': leg, 'type': ev.get('type'),
        't_detect': ev.get('t_detect'), 't_star': ev.get('t_star'),
        'sigma': ev.get('sigma'), 'R2': ev.get('R2'), 'C': ev.get('C'),
        'z_star_re': zs[0], 'z_star_im': zs[1],
        'n_fit_pts': ev.get('n_fit_pts'),
        'labels': ';'.join(ev.get('labels', [])),
        'G_a': g.get('G_a', {}).get('pass'),
        'G_b': g.get('G_b', {}).get('pass'),
        'G_c': g.get('G_c', {}).get('pass'),
        'verdict': ev.get('verdict'),
        'artifact_reason': ev.get('artifact_reason', ''),
    }


def write_figure_data(outdir, cases):
    """Figure data JSONs (prereg §11). cases: {name: case_doc}."""
    os.makedirs(outdir, exist_ok=True)
    # Fig 1: zeta event map (t vs z, type-encoded)
    zc = cases.get('zeta')
    if zc:
        E = validated_events(zc)
        fig1 = [{'type': e['type'], 't_star': e['t_star'],
                 'z_re': e['z_star'][0], 'z_im': e['z_star'][1],
                 'sigma': e['sigma']} for e in E]
        json.dump(fig1, open(os.path.join(outdir, 'fig1_event_map.json'), 'w'))
    # Fig 2: fold/landing fits with residuals
    fits = []
    for name, doc in cases.items():
        for ev in validated_events(doc):
            fits.append({'case': name, 'type': ev['type'],
                         't_star': ev['t_star'], 'sigma': ev.get('sigma'),
                         'R2': ev.get('R2'), 'C': ev.get('C'),
                         'fit_pts': ev.get('fit_pts')})
    json.dump(fits, open(os.path.join(outdir, 'fig2_fits.json'), 'w'))
    # Fig 3: N_off(t) for all cases (down legs)
    noff = {}
    for name, doc in cases.items():
        for lname, leg in doc['legs'].items():
            if 'down' in lname:
                noff[name] = n_off_curve(leg)
    json.dump(noff, open(os.path.join(outdir, 'fig3_noff.json'), 'w'))
    # Fig 4: artifact gallery
    gal = []
    for name, doc in cases.items():
        for ev in all_triggers(doc):
            if ev.get('verdict') in ('artifact', 'suspect'):
                gal.append(event_row(name, '', ev))
    json.dump(gal, open(os.path.join(outdir, 'fig4_artifacts.json'), 'w'))
    return ['fig1_event_map.json', 'fig2_fits.json',
            'fig3_noff.json', 'fig4_artifacts.json']
