#!/usr/bin/env python3
"""
track.py — wide-window Newton tracker for H_t zeros (preregistration v2.0, §§5-6).

Covers: t=0 grid-seeded zero finding (|z| <= 300), bidirectional t-continuation
with adaptive step control, event detection (Type-R / Type-C / Type-X) with
fold/landing fits, and the per-event genuineness protocol (G-a / G-b / G-c).

All tolerances, thresholds and rules are taken literally from prereg §§5-6:
  - Newton: tol |dz| < 1e-30, max 25 iters, require |H| < 1e-25 at convergence.
  - Step: dt0 = 0.01; halve if Newton iters > 12 (floor 0.001);
    double if iters < 4 (cap 0.02); reject |dz| > 0.5 * neighboring gap,
    halve dt and retry.
  - Trigger R: real gap g_k < 0.03; fit g^2 = C^2 (t - t*) on
    g in [1e-7, 3e-2], >= 5 pts, R^2 >= 0.999.
  - Trigger C: complex-pair |Im| < 0.03; linear im^2 fit, same quality rules.
  - Trigger X: non-conjugate off-axis approach < 0.03; logged, never fitted.
  - G-a: >= 10 smooth steps each side of the event (smooth = per-step |dz|
    within the §5 rejection bound).
  - G-b: recompute at J=512; require |dt*| < 3*sigma(fit), |dz*| < 0.05.
  - G-c: |dH/dt + d2H/dz2| < 1e-20 at (z*, t*).

Trajectory JSONs are MANDATORY (§5): every step, every zero records
(t, z, |H|, newton_iters, dt). A case without saved trajectories is VOID.
"""
import json
import math
from mpmath import mp

# ---- preregistered constants (§5) ----
NEWTON_TOL = mp.mpf('1e-30')
NEWTON_MAXIT = 25
H_TOL = mp.mpf('1e-25')
DT0 = mp.mpf('0.01')
DT_MIN = mp.mpf('0.001')
DT_MAX = mp.mpf('0.02')
G_TRIGGER = mp.mpf('0.03')
G_FIT_LO = mp.mpf('1e-7')
G_FIT_HI = mp.mpf('0.03')
IM_TRIGGER = mp.mpf('0.03')
X_TRIGGER = mp.mpf('0.03')
REAL_TOL = mp.mpf('1e-9')      # |Im z| below this => treated as real
OFF_TOL = mp.mpf('1e-6')       # |Im z| above this => counts for N_off
FIT_MIN_PTS = 5
FIT_MIN_R2 = 0.999
DEDUP_DIST = mp.mpf('1e-3')    # pairwise dedup distance for t=0 seeding
GEN_A_STEPS = 10              # G-a: smooth steps required each side
GEN_B_DT_SIGMAS = 3.0
GEN_B_DZ = mp.mpf('0.05')
GEN_C_PDE = mp.mpf('1e-20')


def newton(flow, z0, t, tol=NEWTON_TOL, maxit=NEWTON_MAXIT):
    """Complex Newton on flow.H(., t). Returns (z, iters) or (None, iters).

    Convergence requires |dz| < tol AND |H(z,t)| < H_TOL (prereg §5).
    """
    z = mp.mpc(z0)
    t = mp.mpf(t)
    for it in range(1, maxit + 1):
        f = flow.H(z, t)
        fp = flow.Hp(z, t)
        if fp == 0:
            return None, it
        dz = f / fp
        z -= dz
        if abs(dz) < tol:
            return (z, it) if abs(flow.H(z, t)) < H_TOL else (None, it)
    return None, maxit


def seed_zeros(flow, zmax=300.0, re_step=2.0,
               im_offsets=(0.0, 0.5, -0.5, 1.5, -1.5, 4.0, -4.0)):
    """Grid-seeded complex Newton at t=0. Returns sorted list of mpc zeros
    with |z| <= zmax, deduplicated by pairwise distance (prereg §5)."""
    t0 = mp.mpf('0')
    found = []

    def add(z):
        for w in found:
            if abs(z - w) < DEDUP_DIST:
                return
        found.append(z)

    re = -zmax
    while re <= zmax:
        for im in im_offsets:
            z, _ = newton(flow, mp.mpc(re, im), t0)
            if z is not None and abs(z) <= zmax:
                add(z)
        re += re_step
    # sort: by |z|, then angle, for deterministic labels
    found.sort(key=lambda z: (abs(z), math.atan2(float(mp.im(z)), float(mp.re(z)))))
    return found


def label_zeros(zeros):
    """Deterministic labels: r+{i}/r-{i} for real (|Im|<=REAL_TOL, i by |Re|
    ascending), c+{i}/c-{i} for complex (i by (|Re|,|Im|) ascending)."""
    reals = sorted([z for z in zeros if abs(mp.im(z)) <= REAL_TOL],
                   key=lambda z: abs(mp.re(z)))
    cx = sorted([z for z in zeros if abs(mp.im(z)) > REAL_TOL],
                key=lambda z: (abs(mp.re(z)), abs(mp.im(z))))
    labels = {}
    rp = [z for z in reals if mp.re(z) >= 0]
    rm = [z for z in reals if mp.re(z) < 0]
    for i, z in enumerate(rp, 1):
        labels['r+%d' % i] = z
    for i, z in enumerate(rm, 1):
        labels['r-%d' % i] = z
    # complex: pair conjugates for stable labels
    used = set()
    i = 0
    for z in cx:
        if id(z) in used:
            continue
        i += 1
        labels['c+%d' % i] = z
        used.add(id(z))
        # find conjugate partner
        for w in cx:
            if id(w) in used:
                continue
            if abs(w - mp.conj(z)) < mp.mpf('1e-6'):
                labels['c-%d' % i] = w
                used.add(id(w))
                break
    return labels


def fit_fold(pts):
    """Linear fit of y = x^2 vs t for pts=[(t, x)]. Returns
    (t_star, sigma, R2, C) or None. Prereg §5 quality rules."""
    pts = [(float(t), float(x)) for t, x in pts
           if G_FIT_LO <= mp.mpf(str(abs(float(x)))) <= G_FIT_HI]
    n = len(pts)
    if n < FIT_MIN_PTS:
        return None
    ts = [t for t, _ in pts]
    ys = [x * x for _, x in pts]
    mt = sum(ts) / n
    my = sum(ys) / n
    sxx = sum((t - mt) ** 2 for t in ts)
    sxy = sum((t - mt) * (y - my) for t, y in zip(ts, ys))
    if sxx == 0 or abs(sxy) < 1e-300:
        return None
    slope = sxy / sxx
    icept = my - slope * mt
    t_star = -icept / slope
    ypred = [slope * t + icept for t in ts]
    ss_res = sum((y - yp) ** 2 for y, yp in zip(ys, ypred))
    ss_tot = sum((y - my) ** 2 for y in ys)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    if R2 < FIT_MIN_R2:
        return None
    dof = n - 2
    if dof <= 0:
        return None
    var = ss_res / dof
    sig2 = var / (slope ** 2) * (1.0 / n + mt ** 2 / sxx)
    return (t_star, math.sqrt(sig2), R2, math.sqrt(abs(slope)), n)


def _neighbor_gap(lab, pos):
    z = pos[lab]
    return min((abs(z - pos[o]) for o in pos if o != lab),
               default=mp.mpf('inf'))


def track_leg(flow, seeds, t_start, t_end, dt0=DT0, case_tag="",
             box_seconds=None):
    """t-continuation of labeled zeros from t_start to t_end (prereg §5).

    Returns {'traj', 'events', 'track_fails', 'n_steps', 't_reached',
    'box_expired'}.
    traj[lab] = list of (t, z, |H|, iters, dt_used) — MANDATORY full record.
    events: R/C/X trigger dicts with fit results and raw fit points
    (fit_pts carry (t, z_a, z_b) for the G-b J=512 recompute).
    box_seconds: prereg §10 wall-clock box; on expiry the leg stops early
    and reports the partial map (a preregistered scope cut, not a failure).
    """
    import time as _time
    _t0 = _time.time()
    dt_sign = 1 if mp.mpf(t_end) > mp.mpf(t_start) else -1
    dt = mp.mpf(dt0) * dt_sign
    active = {lab: mp.mpc(z) for lab, z in seeds.items()}
    traj = {lab: [] for lab in active}
    for lab, z in active.items():
        traj[lab].append((mp.mpf(t_start), z, abs(flow.H(z, mp.mpf(t_start))),
                          0, mp.mpf('0')))
    events = []
    track_fails = []
    box_expired = False
    # per-pair accumulation state for R and C triggers
    r_acc = {}   # (la, lb) -> {'pts': [(t, g, za, zb)], 'active': bool}
    c_acc = {}   # frozenset({la, lb}) -> {'pts': [(t, im, z)], 'active': bool}
    x_seen = set()
    t = mp.mpf(t_start)
    n_steps = 0

    def cur_pos():
        return {lab: traj[lab][-1][1] for lab in active}

    while (dt_sign > 0 and t < t_end) or (dt_sign < 0 and t > t_end):
        if box_seconds is not None and _time.time() - _t0 > box_seconds:
            box_expired = True
            break
        t_next = t + dt
        if (dt_sign > 0 and t_next > t_end) or (dt_sign < 0 and t_next < t_end):
            t_next = mp.mpf(t_end)
        pos = cur_pos()
        # Newton step for every active label
        new_pos, iters_map, ok = {}, {}, True
        worst_iters = 0
        for lab, z in active.items():
            nz, it = newton(flow, z, t_next)
            iters_map[lab] = it
            worst_iters = max(worst_iters, it)
            if nz is None:
                ok = False
                break
            if abs(nz - z) > mp.mpf('0.5') * _neighbor_gap(lab, pos):
                ok = False
                break
            new_pos[lab] = nz
        if not ok:
            dt = dt / 2
            if abs(dt) < DT_MIN:
                # prereg §5: halve-and-retry exhausted; drop failed labels,
                # log TRACK_FAIL, continue with the rest (logged, not silent)
                bad = [lab for lab in active if lab not in new_pos]
                for lab in bad:
                    track_fails.append({'t': float(t), 'label': lab,
                                        'case': case_tag})
                    del active[lab]
                if not active:
                    break
                dt = mp.mpf(dt0) * dt_sign
                continue
            continue
        # accept step
        t = t_next
        n_steps += 1
        dt_used = dt
        # adaptive dt for NEXT step (prereg §5)
        if worst_iters > 12:
            dt = max(abs(dt) / 2, DT_MIN) * dt_sign
        elif worst_iters < 4:
            dt = min(abs(dt) * 2, DT_MAX) * dt_sign
        for lab, z in new_pos.items():
            active[lab] = z
            traj[lab].append((t, z, abs(flow.H(z, t)), iters_map[lab],
                              abs(dt_used)))
        if len(active) < 2:
            break
        _detect(flow, t, new_pos, traj, r_acc, c_acc, x_seen, events,
                case_tag, dt_sign)

    return {'traj': traj, 'events': events, 'track_fails': track_fails,
            'n_steps': n_steps, 't_reached': float(t),
            'box_expired': box_expired}


def _detect(flow, t, new_pos, traj, r_acc, c_acc, x_seen, events,
            case_tag, dt_sign):
    """Event detection for one accepted step (prereg §5)."""
    labs = list(new_pos.keys())
    # --- Type-R: adjacent real gaps ---
    real_labs = sorted(
        [lab for lab in labs if abs(mp.im(new_pos[lab])) <= REAL_TOL],
        key=lambda lab: float(mp.re(new_pos[lab])))
    for i in range(len(real_labs) - 1):
        la, lb = real_labs[i], real_labs[i + 1]
        g = abs(mp.re(new_pos[lb]) - mp.re(new_pos[la]))
        key = (la, lb)
        if g < G_TRIGGER:
            acc = r_acc.setdefault(key, {'pts': []})
            acc['pts'].append((t, g, new_pos[la], new_pos[lb]))
        elif key in r_acc and r_acc[key]['pts']:
            # gap reopened without collision: close accumulator, no event
            del r_acc[key]
    # finalize R candidates: gap closed below fit floor, or pair broke
    # (Newton failure -> fit on accumulated points, prereg §5)
    for key in list(r_acc.keys()):
        la, lb = key
        dead = la not in new_pos or lb not in new_pos
        g_now = (abs(mp.re(new_pos[lb]) - mp.re(new_pos[la]))
                 if not dead else None)
        if (g_now is not None and g_now < G_FIT_LO) or dead:
            _finalize_R(r_acc, key, events, case_tag, t)
    # --- Type-C: conjugate pairs' |Im| ---
    pairs = _conjugate_pairs(new_pos)
    for key, (la, lb) in pairs.items():
        im = abs(mp.im(new_pos[la]))
        if im < IM_TRIGGER:
            acc = c_acc.setdefault(key, {'pts': []})
            # store the member with im >= 0
            z = new_pos[la] if mp.im(new_pos[la]) >= 0 else new_pos[lb]
            acc['pts'].append((t, im, z))
        elif key in c_acc and c_acc[key]['pts']:
            del c_acc[key]
    for key in list(c_acc.keys()):
        la, lb = tuple(key)
        dead = la not in new_pos or lb not in new_pos
        im_now = abs(mp.im(new_pos[la])) if not dead else None
        if (im_now is not None and im_now < G_FIT_LO) or dead:
            _finalize_C(c_acc, key, events, case_tag, t)
    # --- Type-X: non-conjugate off-axis approaches ---
    conj_keys = set(pairs.keys())
    for i in range(len(labs)):
        for j in range(i + 1, len(labs)):
            la, lb = labs[i], labs[j]
            if frozenset((la, lb)) in conj_keys:
                continue
            za, zb = new_pos[la], new_pos[lb]
            if abs(mp.im(za)) <= REAL_TOL and abs(mp.im(zb)) <= REAL_TOL:
                continue  # real-real handled by gap monitor
            d = abs(za - zb)
            key = tuple(sorted((la, lb)))
            if d < X_TRIGGER and key not in x_seen:
                x_seen.add(key)
                events.append({'type': 'X', 't_detect': float(t),
                               'labels': [la, lb], 'case': case_tag,
                               'z_a': [float(mp.re(za)), float(mp.im(za))],
                               'z_b': [float(mp.re(zb)), float(mp.im(zb))],
                               'dist': float(d),
                               'note': 'off-axis non-conjugate approach; '
                                       'never fitted (prereg §5)'})


def _conjugate_pairs(pos):
    """Map frozenset({la, lb}) -> (la, lb) for conjugate pairs."""
    out = {}
    labs = list(pos.keys())
    used = set()
    for la in labs:
        if la in used:
            continue
        za = pos[la]
        if abs(mp.im(za)) <= REAL_TOL:
            continue
        for lb in labs:
            if lb == la or lb in used:
                continue
            if abs(pos[lb] - mp.conj(za)) < mp.mpf('1e-6'):
                out[frozenset((la, lb))] = (la, lb)
                used.add(la)
                used.add(lb)
                break
    return out


def _finalize_R(r_acc, key, events, case_tag, t):
    la, lb = key
    pts = r_acc[key]['pts']
    fit = fit_fold([(t_, g_) for t_, g_, _, _ in pts])
    ev = {'type': 'R', 't_detect': float(t), 'labels': [la, lb],
          'case': case_tag, 'n_fit_pts': len(pts),
          'fit_pts': [[float(t_), float(g_),
                       [float(mp.re(za)), float(mp.im(za))],
                       [float(mp.re(zb)), float(mp.im(zb))]]
                      for t_, g_, za, zb in pts]}
    if fit:
        t_star, sig, R2, C, n = fit
        ev.update({'t_star': t_star, 'sigma': sig, 'R2': R2, 'C': C,
                   'z_star': [float((mp.re(pts[-1][2]) + mp.re(pts[-1][3])) / 2),
                              float((mp.im(pts[-1][2]) + mp.im(pts[-1][3])) / 2)]})
    else:
        ev.update({'t_star': None, 'note': 'fit quality rules not met '
                   '(<5 pts or R2<0.999); trigger logged, not mapped'})
    events.append(ev)
    del r_acc[key]


def _finalize_C(c_acc, key, events, case_tag, t):
    la, lb = tuple(key)
    pts = c_acc[key]['pts']
    fit = fit_fold([(t_, im_) for t_, im_, _ in pts])
    ev = {'type': 'C', 't_detect': float(t), 'labels': [la, lb],
          'case': case_tag, 'n_fit_pts': len(pts),
          'fit_pts': [[float(t_), float(im_),
                       [float(mp.re(z)), float(mp.im(z))]]
                      for t_, im_, z in pts]}
    if fit:
        t_star, sig, R2, C, n = fit
        zc = pts[-1][2]
        ev.update({'t_star': t_star, 'sigma': sig, 'R2': R2,
                   'z_star': [float(mp.re(zc)), 0.0]})
    else:
        ev.update({'t_star': None, 'note': 'fit quality rules not met; '
                   'trigger logged, not mapped'})
    events.append(ev)
    del c_acc[key]


# ---- genuineness protocol (prereg §6) ----

def validate_event(flow256, flow512, event, traj):
    """Apply G-a / G-b / G-c to one event dict (in place). Sets
    event['genuineness'] = {'G_a': {...}, 'G_b': {...}, 'G_c': {...}} and
    event['verdict'] in {'genuine', 'suspect', 'artifact'} with
    event['artifact_reason'] when applicable."""
    gen = {}
    # G-a: t-continuation smoothness, >=10 steps each side
    ok_a, det_a = _gen_a(event, traj)
    gen['G_a'] = {'pass': ok_a, 'detail': det_a}
    # G-b: J=512 recompute
    ok_b, det_b = _gen_b(flow512, event)
    gen['G_b'] = {'pass': ok_b, 'detail': det_b}
    # G-c: PDE residual at (z*, t*)
    ok_c, det_c = _gen_c(flow256, event)
    gen['G_c'] = {'pass': ok_c, 'detail': det_c}
    event['genuineness'] = gen
    if event['type'] == 'X':
        event['verdict'] = 'logged'  # X events are never fitted/mapped
        return event
    if event.get('t_star') is None:
        event['verdict'] = 'artifact'
        event['artifact_reason'] = 'fit quality rules not met at detection'
        return event
    if not ok_a:
        event['verdict'] = 'suspect'
        event['artifact_reason'] = 'G-a: ' + det_a
    elif not (ok_b and ok_c):
        event['verdict'] = 'artifact'
        event['artifact_reason'] = 'G-b: %s; G-c: %s' % (det_b, det_c)
    else:
        event['verdict'] = 'genuine'
    return event


def _gen_a(event, traj):
    """>=10 smooth steps each side of t_detect for every participating label.
    Smooth = per-step |dz| within the §5 rejection bound
    (0.5 * distance to nearest other zero at that step)."""
    t_det = event['t_detect']
    for lab in event['labels']:
        pts = traj.get(lab)
        if not pts:
            return False, 'no trajectory for %s' % lab
        # index of step nearest t_detect
        idx = min(range(len(pts)), key=lambda i: abs(pts[i][0] - t_det))
        if idx < GEN_A_STEPS or idx + GEN_A_STEPS >= len(pts):
            return False, ('fewer than %d steps on one side for %s '
                           '(SUSPECT per §6)' % (GEN_A_STEPS, lab))
        for i in range(idx - GEN_A_STEPS, idx + GEN_A_STEPS):
            t_i, z_i = pts[i][0], pts[i][1]
            t_j, z_j = pts[i + 1][0], pts[i + 1][1]
            dz = abs(z_j - z_i)
            gap = _min_dist_to_others(lab, traj, i)
            if gap is None:
                continue
            if dz > mp.mpf('0.5') * gap:
                return False, ('non-smooth step for %s at t=%.4f '
                               '(|dz|=%.2e > 0.5*gap)' % (lab, float(t_i),
                                                          float(dz)))
    return True, 'smooth on [%d, %d] steps each side' % (GEN_A_STEPS,
                                                        GEN_A_STEPS)


def pts_at(pts, t):
    for p in pts:
        if abs(p[0] - t) < mp.mpf('1e-12'):
            return p
    return None


def _min_dist_to_others(lab, traj, i):
    z = traj[lab][i][1]
    best = None
    for l2, pts in traj.items():
        if l2 == lab or i >= len(pts):
            continue
        d = abs(z - pts[i][1])
        if best is None or d < best:
            best = d
    return best


def _gen_b(flow512, event):
    """Recompute the event at J=512: re-polish fit-point zeros, refit,
    require |dt*| < 3*sigma and |dz*| < 0.05 (prereg §6)."""
    if event.get('t_star') is None or not event.get('fit_pts'):
        return False, 'no fit to recompute'
    sig = event.get('sigma')
    if sig is None:
        return False, 'no sigma from fit'
    new_pts = []
    z_a_new = z_b_new = None
    for fp in event['fit_pts']:
        t_ = fp[0]
        t = mp.mpf(str(t_))
        if event['type'] == 'R':
            x_, za_, zb_ = fp[1], fp[2], fp[3]
            za = mp.mpc(za_[0], za_[1])
            zb = mp.mpc(zb_[0], zb_[1])
            na, _ = newton(flow512, za, t)
            nb, _ = newton(flow512, zb, t)
            if na is None or nb is None:
                return False, 'J=512 Newton failed at t=%.4f' % t_
            g = abs(na - nb)
            new_pts.append((t_, float(g)))
            z_a_new, z_b_new = na, nb
        else:  # C: fit_pts entries are [t, im, [re, im]]
            za_ = fp[2]
            z = mp.mpc(za_[0], za_[1])
            nz, _ = newton(flow512, z, t)
            if nz is None:
                return False, 'J=512 Newton failed at t=%.4f' % t_
            new_pts.append((t_, float(abs(mp.im(nz)))))
            z_a_new = nz
    fit = fit_fold(new_pts)
    if not fit:
        return False, 'J=512 refit failed quality rules'
    t_s2, sig2, R2_2, C2, n2 = fit
    dt_star = abs(t_s2 - event['t_star'])
    if dt_star >= GEN_B_DT_SIGMAS * sig:
        return False, '|dt*|=%.2e >= 3*sigma' % dt_star
    if event['type'] == 'R':
        z_new = (z_a_new + z_b_new) / 2
    else:
        z_new = mp.mpc(mp.re(z_a_new), 0)
    z_old = mp.mpc(event['z_star'][0], event['z_star'][1])
    dz_star = abs(z_new - z_old)
    if dz_star >= GEN_B_DZ:
        return False, '|dz*|=%.2e >= 0.05' % float(dz_star)
    return True, 'J=512: dt*=%.2e (<3sig), dz*=%.2e (<0.05), R2=%.5f' % (
        dt_star, float(dz_star), R2_2)


def _gen_c(flow256, event):
    """PDE residual < 1e-20 at (z*, t*) (prereg §6).

    v2.1: uses pde_residual_exact (analytic dH/dt); the v2.0 finite-
    difference version measured its own truncation, not the PDE.
    """
    if event.get('t_star') is None or not event.get('z_star'):
        return False, 'no (z*, t*) to check'
    z = mp.mpc(event['z_star'][0], event['z_star'][1])
    t = mp.mpf(str(event['t_star']))
    pres = getattr(flow256, 'pde_residual_exact', None) or flow256.pde_residual
    try:
        r = pres(z, t)
    except AttributeError:
        return False, 'flow lacks pde_residual'
    ok = r < GEN_C_PDE
    return ok, '|dH/dt + d2H/dz2| = %.2e %s 1e-20' % (
        float(r), '<' if ok else '>=')


# ---- serialization (prereg §11: trajectory JSONs are MANDATORY) ----

def _f(x):
    return float(x)


def traj_to_jsonable(traj):
    return {lab: [[_f(t), _f(mp.re(z)), _f(mp.im(z)), _f(ah), it, _f(dt)]
                  for t, z, ah, it, dt in pts]
            for lab, pts in traj.items()}


def event_to_jsonable(ev):
    out = {}
    for k, v in ev.items():
        if k in ('fit_pts', 'genuineness'):
            out[k] = v
        elif isinstance(v, (mp.mpf, mp.mpc)):
            out[k] = complex(v) if isinstance(v, mp.mpc) else float(v)
        elif isinstance(v, complex):
            out[k] = [v.real, v.imag]
        else:
            out[k] = v
    return out


def save_case_json(path, case, prereg_tag, flow_info, seeds, legs):
    """Write the prereg §11 case file. legs: {name: track_leg result}."""
    doc = {
        'case': case,
        'prereg': prereg_tag,
        'flow': flow_info,
        'seeds': {k: [_f(mp.re(v)), _f(mp.im(v))] for k, v in seeds.items()},
        'legs': {},
    }
    for name, r in legs.items():
        doc['legs'][name] = {
            'traj': traj_to_jsonable(r['traj']),
            'events': [event_to_jsonable(e) for e in r['events']],
            'track_fails': r['track_fails'],
            'n_steps': r['n_steps'],
            't_reached': r['t_reached'],
            'box_expired': r.get('box_expired', False),
        }
    with open(path, 'w') as f:
        json.dump(doc, f)
    return path
